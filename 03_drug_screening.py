#!/usr/bin/env python3
"""
03_drug_screening.py — LINCS L1000 signature-reversal screening for spaceflight
transcriptomic signatures.

Queries the SigCom LINCS API for compounds that REVERSE the consensus spaceflight
signatures (liver/immune and brain), aggregates to per-compound tiers, ranks by
a reproducibility-weighted composite score, and annotates with Broad Drug
Repurposing Hub metadata (MoA, target, clinical phase, SMILES).

Falls back to local GMT enrichment if the SigCom API is unreachable.

Usage:
  python 03_drug_screening.py
"""
import json, os, re, sys, time
import numpy as np
import pandas as pd
import requests

# ── Paths ──────────────────────────────────────────────────────────────────
PROJECT_ROOT = "/workspace/astronaut-opposite-forcing"
PROC_DIR     = os.path.join(PROJECT_ROOT, "data", "processed")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "tables")
FIGURES_DIR  = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── SigCom API ─────────────────────────────────────────────────────────────
METADATA_API = "https://maayanlab.cloud/sigcom-lincs/metadata-api"
DATA_API     = "https://maayanlab.cloud/sigcom-lincs/data-api/api/v1"

# ── Data lake resources ────────────────────────────────────────────────────
GMT_DRUG = ("/mnt/datalake/LINCS1000/RNAseq_transcriptomics_genesets/"
            "single_drug_perturbations-v1.0.gmt")
HUB_MOL  = ("/mnt/datalake/broad_drug_repurposing_hub/"
            "broad_repurposing_hub_molecule_with_smiles.parquet")
HUB_INFO = ("/mnt/datalake/broad_drug_repurposing_hub/"
            "broad_repurposing_hub_phase_moa_target_info.parquet")

# ── Tissue-relevant L1000 cell lines ───────────────────────────────────────
# Liver-relevant: hepatocyte-like lines + immune lines
LIVER_TISSUE_LINES = ["HEPG2", "HEPG2", "HEK293T", "A375", "HA1E", "MCF7",
                      "PC3", "VCAP", "A549", "HT29", "HCT116", "JURKAT",
                      "THP1", "U937", "RAJI", "BT474", "T47D", "MDAMB231"]
# Brain-relevant: neuronal/glial lines
BRAIN_TISSUE_LINES = ["SHSY5Y", "NCCIT", "SKNSH", "YAPC", "NEURO2A",
                      "U87", "U251", "T98G", "LN229", "A172", "C6",
                      "HEK293T", "MCF7", "A375", "HA1E", "PC3", "A549"]

# ── Known spaceflight countermeasure drugs for positive-control recovery ───
# These are compounds with published evidence in spaceflight/microgravity/
# radiation countermeasure contexts
KNOWN_COUNTERMEASURES = [
    "resveratrol", "metformin", "rapamycin", "dexamethasone",
    "NAC", "N-acetylcysteine", "ascorbic acid", "vitamin C",
    "trolox", "tempol", "ebselen", "sulforaphane",
    "curcumin", "quercetin", "ibuprofen", "aspirin",
    "lithium", "valproic acid", "trichostatin A", "SAHA",
    "vorinostat", "entinostat", "panobinostat",
]


# ═══════════════════════════════════════════════════════════════════════════
# IO helpers
# ═══════════════════════════════════════════════════════════════════════════
def read_genes(path):
    """Read gene symbols from a text file (one per line)."""
    with open(path) as f:
        toks = re.split(r"[\s,;]+", f.read().strip())
    return sorted({t.strip().upper() for t in toks if t.strip()})


# ═══════════════════════════════════════════════════════════════════════════
# SigCom API functions
# ═══════════════════════════════════════════════════════════════════════════
def api_reachable(timeout=15):
    """Probe with the real gene-resolution endpoint."""
    try:
        r = requests.post(
            f"{METADATA_API}/entities/find",
            json={"filter": {"where": {"meta.symbol": {"inq": ["TNF"]}}}},
            timeout=timeout
        )
        return r.status_code == 200
    except Exception:
        return False


def resolve_genes(genes, timeout=60):
    """Resolve HGNC symbols to L1000 entity UUIDs."""
    r = requests.post(
        f"{METADATA_API}/entities/find",
        json={"filter": {"where": {"meta.symbol": {"inq": list(genes)}}}},
        timeout=timeout
    )
    r.raise_for_status()
    return {e["meta"]["symbol"]: e["id"] for e in r.json()}


def connectivity_query(up_ids, dn_ids, database="l1000_cp", limit=2000, timeout=300,
                       max_retries=5, initial_wait=10):
    """Two-sided reversal query against L1000, with retry on 503."""
    body = {
        "up_entities": list(up_ids),
        "down_entities": list(dn_ids),
        "limit": limit,
        "database": database
    }
    for attempt in range(max_retries):
        try:
            r = requests.post(f"{DATA_API}/enrich/ranktwosided", json=body, timeout=timeout)
            if r.status_code == 200:
                return r.json()
            elif r.status_code == 503:
                wait = initial_wait * (attempt + 1)
                print(f"    503 (no available server) — retry {attempt+1}/{max_retries} in {wait}s...")
                time.sleep(wait)
            else:
                r.raise_for_status()
        except requests.exceptions.Timeout:
            wait = initial_wait * (attempt + 1)
            print(f"    Timeout — retry {attempt+1}/{max_retries} in {wait}s...")
            time.sleep(wait)
    # All retries exhausted
    print(f"    All {max_retries} retries exhausted. API unavailable.")
    return None


def signature_meta(uuids, batch=100, timeout=120):
    """Resolve signature UUIDs to compound metadata."""
    out = {}
    for i in range(0, len(uuids), batch):
        chunk = uuids[i:i+batch]
        r = requests.post(
            f"{METADATA_API}/signatures/find",
            json={"filter": {"where": {"id": {"inq": chunk}}}},
            timeout=timeout
        )
        r.raise_for_status()
        for s in r.json():
            m = s.get("meta", {})
            out[s["id"]] = {k: m.get(k) for k in
                            ("pert_name", "cell_line", "pert_dose", "pert_time",
                             "pubchem_id", "cmap_id", "moa")}
    return out


# ═══════════════════════════════════════════════════════════════════════════
# BRD → compound name mapping (Broad Drug Repurposing Hub)
# ═══════════════════════════════════════════════════════════════════════════
def base_brd(bid):
    m = re.match(r"(BRD-[A-Z0-9]+)", str(bid))
    return m.group(1) if m else None


def load_brd_map():
    """Load BRD → pert_iname mapping from Broad DRH."""
    if not os.path.exists(HUB_MOL):
        print(f"  WARNING: Broad DRH molecule file not found at {HUB_MOL}")
        return {}, None
    hub = pd.read_parquet(HUB_MOL)
    brd_col = next((c for c in hub.columns
                    if re.search(r"broad|brd|deprecated_broad_id|pert_id", c, re.I)), None)
    name_col = next((c for c in hub.columns
                     if re.search(r"pert_iname|name", c, re.I)), None)
    m = {}
    if brd_col and name_col:
        for b, n in zip(hub[brd_col], hub[name_col]):
            bb = base_brd(b)
            if bb and isinstance(n, str) and n:
                m.setdefault(bb, n)
    print(f"  Loaded {len(m)} BRD → name mappings from Broad DRH")
    return m, hub


def load_hub_info():
    """Load MoA/target/phase info from Broad DRH."""
    if not os.path.exists(HUB_INFO):
        return None
    info = pd.read_parquet(HUB_INFO)
    print(f"  Loaded Broad DRH phase/MoA/target info: {info.shape}")
    return info


def annotate_name(pert_name, brd_map):
    if isinstance(pert_name, str) and pert_name and not pert_name.startswith("BRD-"):
        return pert_name
    bb = base_brd(pert_name)
    return brd_map.get(bb, pert_name)


# ═══════════════════════════════════════════════════════════════════════════
# Aggregation + scoring
# ═══════════════════════════════════════════════════════════════════════════
def zscore(s):
    s = np.asarray(s, float)
    return (s - np.nanmean(s)) / (np.nanstd(s) + 1e-9)


def build_reverser_table(results, meta, brd_map):
    df = pd.DataFrame(results)
    df["uuid"] = df.get("uuid", df.get("id"))
    df["compound_raw"] = df["uuid"].map(lambda u: (meta.get(u, {}) or {}).get("pert_name"))
    df["cell_line"]    = df["uuid"].map(lambda u: (meta.get(u, {}) or {}).get("cell_line"))
    df["moa"]          = df["uuid"].map(lambda u: (meta.get(u, {}) or {}).get("moa"))
    df["compound"]     = df["compound_raw"].map(lambda x: annotate_name(x, brd_map))
    return df


def aggregate(df, min_sigs=2):
    rev = df[df["type"] == "reversers"].copy()
    mim = df[df["type"] == "mimickers"].copy()
    n_mim = mim.groupby("compound").size().rename("n_mimicking_sigs")
    g = rev.groupby("compound")
    agg = pd.DataFrame({
        "n_reversing_sigs": g.size(),
        "n_cell_lines":     g["cell_line"].nunique(),
        "median_z_sum":     g["z-sum"].median(),
        "best_z_sum":       g["z-sum"].min(),
        "best_fdr_down":    g["fdr-down"].min(),
        "moa":              g["moa"].agg(lambda s: next((x for x in s if isinstance(x, str) and x), None)),
    }).join(n_mim).fillna({"n_mimicking_sigs": 0})
    agg["reverser_specificity"] = agg["n_reversing_sigs"] / (agg["n_reversing_sigs"] + agg["n_mimicking_sigs"])
    agg = agg.reset_index()

    tier1 = agg[agg["n_reversing_sigs"] >= min_sigs].copy()
    tier2 = agg[agg["n_reversing_sigs"] < min_sigs].copy().sort_values("best_z_sum")

    if not tier1.empty:
        strength = -tier1["median_z_sum"]
        repro    = np.log1p(tier1["n_reversing_sigs"]) + np.log1p(tier1["n_cell_lines"])
        signif   = -np.log10(tier1["best_fdr_down"].clip(lower=1e-320))
        tier1["reverser_score"] = (0.45 * zscore(repro) + 0.30 * zscore(strength) + 0.25 * zscore(signif)) \
                                  * tier1["reverser_specificity"]
        tier1 = tier1.sort_values("reverser_score", ascending=False).reset_index(drop=True)
        tier1.insert(0, "rank", tier1.index + 1)

    return tier1, tier2, rev, mim


# ═══════════════════════════════════════════════════════════════════════════
# Positive-control recovery
# ═══════════════════════════════════════════════════════════════════════════
def check_positive_controls(tier1, tier2, known_drugs):
    """Check which known countermeasure drugs appear in the reverser rankings."""
    results = []
    all_compounds_t1 = set(tier1["compound"].dropna().str.lower()) if not tier1.empty else set()
    all_compounds_t2 = set(tier2["compound"].dropna().str.lower()) if not tier2.empty else set()

    for drug in known_drugs:
        drug_lower = drug.lower()
        in_t1 = drug_lower in all_compounds_t1
        in_t2 = drug_lower in all_compounds_t2
        rank_val = None
        if in_t1:
            match = tier1[tier1["compound"].str.lower() == drug_lower]
            if not match.empty:
                rank_val = int(match["rank"].iloc[0])
        results.append({
            "drug": drug,
            "in_tier1": in_t1,
            "in_tier2": in_t2,
            "tier1_rank": rank_val,
            "recovered": in_t1 or in_t2
        })
    return pd.DataFrame(results)


# ═══════════════════════════════════════════════════════════════════════════
# Local GMT fallback
# ═══════════════════════════════════════════════════════════════════════════
def parse_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sets[parts[0]] = {g.split(",")[0].strip().upper() for g in parts[2:] if g.strip()}
    return sets


def fallback_local(up, dn, out_dir, min_sigs=2, universe=12000):
    from scipy.stats import hypergeom
    if not os.path.exists(GMT_DRUG):
        raise SystemExit(f"API unreachable and local GMT not found at {GMT_DRUG}")
    sets = parse_gmt(GMT_DRUG)
    up, dn = set(up), set(dn)
    rows = []
    for name, genes in sets.items():
        low = name.lower()
        is_up = ("-up" in low) or low.endswith(" up") or ("_up" in low)
        is_dn = ("-dn" in low) or ("down" in low) or ("_dn" in low)
        if not (is_up or is_dn):
            continue
        q = dn if is_up else up
        k = len(genes & q)
        if k == 0:
            continue
        p = hypergeom.sf(k - 1, universe, len(q), len(genes))
        rows.append({"signature": name, "compound": re.split(r"[-_ ]", name)[0],
                     "overlap": k, "p": p, "direction": "up" if is_up else "dn"})
    fb = pd.DataFrame(rows)
    if fb.empty:
        raise SystemExit("Local fallback produced no overlaps; check gene symbols.")
    fb["z-sum"] = np.log10(fb["p"].clip(lower=1e-320))
    g = fb.groupby("compound")
    agg = pd.DataFrame({
        "n_reversing_sigs": g.size(),
        "median_z_sum": g["z-sum"].median(),
        "best_z_sum": g["z-sum"].min(),
        "best_p": g["p"].min()
    }).reset_index()
    agg = agg.sort_values("best_z_sum").reset_index(drop=True)
    agg.insert(0, "rank", agg.index + 1)
    agg["method"] = "local_gmt_fallback"
    fb.to_csv(os.path.join(out_dir, "reversers_raw.csv"), index=False)
    agg.to_csv(os.path.join(out_dir, "tier1_ranking.csv"), index=False)
    return agg, fb, "local_gmt_fallback"


# ═══════════════════════════════════════════════════════════════════════════
# Main screening function
# ═══════════════════════════════════════════════════════════════════════════
def screen_signature(stratum, up_file, down_file, tissue_lines, out_dir):
    """Run full LINCS screening pipeline for one tissue signature."""
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{'='*70}")
    print(f"  Screening: {stratum}")
    print(f"{'='*70}")

    up_genes = read_genes(up_file)
    dn_genes = read_genes(down_file)
    print(f"  Signature: {len(up_genes)} up, {len(dn_genes)} down genes")

    if len(up_genes) == 0 or len(dn_genes) == 0:
        print(f"  ERROR: Empty gene list for {stratum}")
        return None

    # Check API reachability
    if not api_reachable():
        print("  SigCom API unreachable → local GMT fallback")
        agg, raw, method = fallback_local(up_genes, dn_genes, out_dir)
        summary = {
            "stratum": stratum, "engine": method,
            "signature_up_genes": len(up_genes), "signature_dn_genes": len(dn_genes),
            "unique_reverser_compounds": int(agg.shape[0]),
            "note": "API unreachable; indicative ranks only"
        }
        with open(os.path.join(out_dir, "robustness_summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        return summary

    # Resolve genes to L1000 entities
    print("  Resolving genes to L1000 entities...")
    sym2id = resolve_genes(up_genes + dn_genes)
    up_ids = [sym2id[g] for g in up_genes if g in sym2id]
    dn_ids = [sym2id[g] for g in dn_genes if g in sym2id]
    cov_up = 100.0 * len(up_ids) / max(len(up_genes), 1)
    cov_dn = 100.0 * len(dn_ids) / max(len(dn_genes), 1)
    print(f"  L1000 coverage: up={cov_up:.1f}% ({len(up_ids)}/{len(up_genes)}), "
          f"down={cov_dn:.1f}% ({len(dn_ids)}/{len(dn_genes)})")

    unresolved_up = [g for g in up_genes if g not in sym2id]
    unresolved_dn = [g for g in dn_genes if g not in sym2id]
    if unresolved_up:
        print(f"  Unresolved up genes (first 10): {unresolved_up[:10]}")
    if unresolved_dn:
        print(f"  Unresolved down genes (first 10): {unresolved_dn[:10]}")

    # Connectivity query
    print("  Querying SigCom LINCS (l1000_cp)...")
    resp = connectivity_query(up_ids, dn_ids, database="l1000_cp", limit=2000,
                              max_retries=5, initial_wait=15)

    if resp is None:
        print("  SigCom data API unavailable → local GMT fallback")
        agg, raw, method = fallback_local(up_genes, dn_genes, out_dir)
        summary = {
            "stratum": stratum, "engine": method,
            "signature_up_genes": len(up_genes), "signature_dn_genes": len(dn_genes),
            "l1000_coverage_up_pct": round(cov_up, 1),
            "l1000_coverage_dn_pct": round(cov_dn, 1),
            "unique_reverser_compounds": int(agg.shape[0]),
            "note": "SigCom data API returned 503; local GMT fallback used"
        }
        with open(os.path.join(out_dir, "robustness_summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        return summary

    n_rev = resp.get("reversers")
    n_mim = resp.get("mimickers")
    print(f"  Results: {len(resp['results'])} rows | reversers={n_rev} mimickers={n_mim}")

    # Resolve signature metadata
    print("  Resolving signature metadata...")
    uuids = [r.get("uuid", r.get("id")) for r in resp["results"]]
    meta = signature_meta(uuids)

    # Build reverser table
    brd_map, hub_mol = load_brd_map()
    df = build_reverser_table(resp["results"], meta, brd_map)

    # Aggregate + tier
    tier1, tier2, rev, mim = aggregate(df, min_sigs=2)
    print(f"  Tier-1 (reproducible): {tier1.shape[0]} compounds")
    print(f"  Tier-2 (single-sig):   {tier2.shape[0]} compounds")

    # Annotate with Broad DRH phase/MoA/target
    hub_info = load_hub_info()
    if hub_info is not None and not tier1.empty:
        # Try to merge on compound name
        name_col = next((c for c in hub_info.columns
                         if re.search(r"pert_iname|name", c, re.I)), None)
        if name_col:
            # Deduplicate hub_info by name (keep first)
            hub_info_dedup = hub_info.drop_duplicates(subset=[name_col])
            tier1 = tier1.merge(
                hub_info_dedup.add_prefix("drh_"),
                left_on="compound", right_on=f"drh_{name_col}",
                how="left"
            )

    # Tissue-context view
    tissue_set = {t.strip().upper() for t in tissue_lines}
    rev["cell_up"] = rev["cell_line"].astype(str).str.upper()
    tissue_hits = sorted(
        set(rev.loc[rev["cell_up"].isin(tissue_set), "compound"].dropna())
        & set(tier1["compound"]) if not tier1.empty else set()
    )
    print(f"  Tissue-supported Tier-1: {len(tissue_hits)} compounds")

    # Positive-control recovery
    pc_results = check_positive_controls(tier1, tier2, KNOWN_COUNTERMEASURES)
    n_recovered = pc_results["recovered"].sum()
    print(f"  Positive-control recovery: {n_recovered}/{len(KNOWN_COUNTERMEASURES)} known drugs found")

    # Save outputs
    df.to_csv(os.path.join(out_dir, "reversers_raw.csv"), index=False)
    tier1.to_csv(os.path.join(out_dir, "tier1_ranking.csv"), index=False)
    tier2.to_csv(os.path.join(out_dir, "tier2_single_signature.csv"), index=False)
    pc_results.to_csv(os.path.join(out_dir, "positive_control_recovery.csv"), index=False)

    # Top 20 reversers
    print(f"\n  Top 20 Tier-1 reversers for {stratum}:")
    if not tier1.empty:
        display_cols = [c for c in ["rank", "compound", "n_reversing_sigs", "n_cell_lines",
                                     "median_z_sum", "reverser_score", "moa"]
                        if c in tier1.columns]
        print(tier1[display_cols].head(20).to_string(index=False))

    # Robustness summary
    summary = {
        "stratum": stratum,
        "engine": "sigcom_lincs",
        "database": "l1000_cp",
        "signature_up_genes": len(up_genes),
        "signature_dn_genes": len(dn_genes),
        "l1000_coverage_up_pct": round(cov_up, 1),
        "l1000_coverage_dn_pct": round(cov_dn, 1),
        "db_reversers_count": int(n_rev) if n_rev is not None else None,
        "db_mimickers_count": int(n_mim) if n_mim is not None else None,
        "rows_retrieved": int(len(resp["results"])),
        "unique_reverser_compounds": int(rev["compound"].nunique()),
        "tier1_reproducible_compounds": int(tier1.shape[0]),
        "tier2_single_signature": int(tier2.shape[0]),
        "strongest_reverser_zsum": float(tier1["median_z_sum"].min()) if not tier1.empty else None,
        "tissue_lines_used": sorted(tissue_set),
        "tier1_tissue_supported": tissue_hits,
        "tier1_tissue_supported_n": len(tissue_hits),
        "positive_control_recovery": {
            "n_known_drugs": len(KNOWN_COUNTERMEASURES),
            "n_recovered": int(n_recovered),
            "recovered_drugs": pc_results[pc_results["recovered"]]["drug"].tolist(),
        }
    }
    with open(os.path.join(out_dir, "robustness_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    return summary


# ═══════════════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("  LINCS L1000 Signature-Reversal Drug Screening")
    print("  Spaceflight Transcriptomic Signatures → Opposite Forcings")
    print("=" * 70)

    # Screen both tissue signatures
    summaries = {}

    # Liver/immune signature
    li_up = os.path.join(PROC_DIR, "signature_liver_immune_up_genes.txt")
    li_dn = os.path.join(PROC_DIR, "signature_liver_immune_down_genes.txt")
    li_out = os.path.join(RESULTS_DIR, "lincs_liver_immune")
    if os.path.exists(li_up) and os.path.exists(li_dn):
        summaries["liver_immune"] = screen_signature(
            "liver_immune", li_up, li_dn, LIVER_TISSUE_LINES, li_out
        )
    else:
        print(f"  WARNING: Liver/immune signature files not found")

    # Brain signature
    br_up = os.path.join(PROC_DIR, "signature_brain_up_genes.txt")
    br_dn = os.path.join(PROC_DIR, "signature_brain_down_genes.txt")
    br_out = os.path.join(RESULTS_DIR, "lincs_brain")
    if os.path.exists(br_up) and os.path.exists(br_dn):
        summaries["brain"] = screen_signature(
            "brain", br_up, br_dn, BRAIN_TISSUE_LINES, br_out
        )
    else:
        print(f"  WARNING: Brain signature files not found")

    # Combined summary
    with open(os.path.join(RESULTS_DIR, "lincs_screening_summary.json"), "w") as f:
        json.dump(summaries, f, indent=2)

    print(f"\n{'='*70}")
    print("  Screening complete!")
    print(f"{'='*70}")
    for stratum, s in summaries.items():
        if s:
            print(f"  {stratum}: {s.get('tier1_reproducible_compounds', 'N/A')} Tier-1 compounds, "
                  f"coverage up={s.get('l1000_coverage_up_pct', 'N/A')}% "
                  f"down={s.get('l1000_coverage_dn_pct', 'N/A')}%")
    print(f"\n  Results saved to: {RESULTS_DIR}/")


if __name__ == "__main__":
    main()
