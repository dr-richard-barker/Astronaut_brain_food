#!/usr/bin/env python3
"""
05_nutrient_gene_mapping.py — Map LINCS drug hits to food compounds and vegan sources.

Creates a nutrient-gene-therapeutic bridge:
1. Load LINCS reverser compounds (from Script 03)
2. Match compounds to FooDB food compounds (by name)
3. Identify vegan food sources for each matched compound
4. Map food compounds to signature genes via LINCS drug targets
"""
import json, os, re, sys
import pandas as pd
import numpy as np

PROJECT_ROOT = "/workspace/astronaut-opposite-forcing"
PROC_DIR     = os.path.join(PROJECT_ROOT, "data", "processed")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "tables")
FOODB_DIR    = os.path.join(PROJECT_ROOT, "data", "raw", "foodb", "foodb_2020_04_07_csv")

VEGAN_KEYWORDS = [
    "broccoli", "spinach", "kale", "cabbage", "carrot", "tomato", "potato",
    "onion", "garlic", "ginger", "pepper", "cucumber", "lettuce", "celery",
    "asparagus", "artichoke", "beet", "radish", "turnip", "cauliflower",
    "brussels", "zucchini", "eggplant", "mushroom", "fennel", "leek",
    "apple", "orange", "banana", "grape", "berry", "strawberry", "blueberry",
    "raspberry", "cranberry", "lemon", "lime", "grapefruit", "pineapple",
    "mango", "papaya", "peach", "pear", "plum", "cherry", "apricot",
    "pomegranate", "watermelon", "cantaloupe", "fig", "date", "kiwi",
    "avocado", "coconut", "olive",
    "bean", "lentil", "chickpea", "pea", "soybean", "tofu", "tempeh",
    "mung", "adzuki", "black bean", "kidney bean", "navy bean",
    "almond", "walnut", "cashew", "pistachio", "pecan", "hazelnut",
    "brazil nut", "macadamia", "peanut", "sunflower seed", "pumpkin seed",
    "sesame", "flax", "chia", "hemp seed", "pine nut",
    "rice", "wheat", "oat", "barley", "quinoa", "buckwheat", "millet",
    "sorghum", "rye", "corn", "maize", "amaranth", "spelt",
    "turmeric", "cumin", "coriander", "curry", "basil", "oregano",
    "thyme", "rosemary", "sage", "mint", "parsley", "cilantro",
    "dill", "cinnamon", "cloves", "nutmeg", "cardamom", "saffron",
    "paprika", "cayenne",
    "olive oil", "coconut oil", "sesame oil", "vinegar", "soy sauce",
    "mustard", "cocoa", "chocolate", "vanilla",
    "seaweed", "spirulina", "chlorella", "nori", "wakame", "kelp",
    "green tea", "black tea", "coffee", "cider",
    "capers", "endive", "arugula", "watercress", "radicchio",
    "shallot", "scallion", "chive", "horseradish", "wasabi",
    "star anise", "fennel", "caraway", "anise", "fenugreek",
    "goji", "acai", "mulberry", "elderberry", "blackberry",
    "cabbage", "bok choy", "water chestnut", "lotus",
    "taro", "yam", "sweet potato", "pumpkin", "squash",
]

NON_VEGAN_KEYWORDS = [
    "milk", "cheese", "butter", "cream", "yogurt", "whey", "casein",
    "beef", "pork", "chicken", "turkey", "duck", "goose", "lamb",
    "veal", "mutton", "bacon", "ham", "sausage", "meat",
    "fish", "salmon", "tuna", "cod", "herring", "mackerel", "sardine",
    "shrimp", "prawn", "crab", "lobster", "oyster", "clam", "mussel",
    "egg", "honey", "gelatin", "lard", "tallow", "stock", "broth",
]

KNOWN_NUTRIENTS = {
    'quercetin', 'sulforaphane', 'curcumin', 'resveratrol',
    'ascorbic acid', 'ascorbic', 'vitamin c', 'vitamin e',
    'luteolin', 'apigenin', 'kaempferol', 'catechin',
    'epigallocatechin', 'genistein', 'daidzein',
    'allicin', 'capsaicin', 'piperine', 'cinnamaldehyde',
    'gingerol', 'zingerone', 'thymol', 'carvacrol',
    'eugenol', 'rosmarinic acid', 'carnosic acid',
    'ellagic acid', 'gallic acid', 'ferulic acid',
    'chlorogenic acid', 'caffeic acid', 'coumaric acid',
    'sulforaphane', 'indole-3-carbinol', 'diindolylmethane',
    'berberine', 'fisetin', 'piceatannol', 'pterostilbene',
    'naringenin', 'hesperidin', 'rutin', 'diosmetin',
    'tangeretin', 'nobiletin', 'silymarin', 'silibinin',
    'baicalein', 'wogonin', 'orsellinic',
    'carnosol', 'ursolic acid', 'oleanolic acid',
    'withaferin', 'withanolide', 'ginsenoside',
    'sulforaphane', 'phenethyl isothiocyanate',
    'lycopene', 'beta-carotene', 'beta carotene',
    'lutein', 'zeaxanthin', 'astaxanthin',
    'tocopherol', 'tocotrienol',
}


def is_vegan_food(food_name):
    if not isinstance(food_name, str):
        return False
    name_lower = food_name.lower()
    for kw in NON_VEGAN_KEYWORDS:
        if kw in name_lower:
            return False
    for kw in VEGAN_KEYWORDS:
        if kw in name_lower:
            return True
    return False


def load_foodb_compounds():
    """Load FooDB Compound.csv."""
    fpath = os.path.join(FOODB_DIR, "Compound.csv")
    if not os.path.exists(fpath):
        print(f"  FooDB Compound.csv not found at {fpath}")
        return None
    print(f"  Loading FooDB Compound.csv...")
    df = pd.read_csv(fpath, low_memory=False)
    print(f"    {df.shape[0]} compounds")
    df['name_lower'] = df['name'].str.lower().str.strip()
    return df


def load_foodb_foods():
    """Load FooDB Food.csv."""
    fpath = os.path.join(FOODB_DIR, "Food.csv")
    if not os.path.exists(fpath):
        print(f"  FooDB Food.csv not found at {fpath}")
        return None
    print(f"  Loading FooDB Food.csv...")
    df = pd.read_csv(fpath, low_memory=False)
    print(f"    {df.shape[0]} foods")
    return df


def load_foodb_content(compound_ids):
    """Load FooDB Content.csv filtered to our compound IDs."""
    fpath = os.path.join(FOODB_DIR, "Content.csv")
    if not os.path.exists(fpath):
        print(f"  FooDB Content.csv not found at {fpath}")
        return None
    print(f"  Loading FooDB Content.csv (filtered to {len(compound_ids)} compounds)...")
    # Read in chunks and filter
    chunks = []
    for chunk in pd.read_csv(fpath, chunksize=200000, low_memory=False):
        filtered = chunk[chunk['source_id'].isin(compound_ids)]
        if not filtered.empty:
            chunks.append(filtered)
    if not chunks:
        print(f"    No content records matched")
        return None
    df = pd.concat(chunks, ignore_index=True)
    print(f"    {df.shape[0]} content records matched")
    return df


def match_compounds_to_foodb(lincs_compounds, foodb_compound):
    """Match LINCS hit compounds to FooDB compounds by name."""
    if foodb_compound is None:
        return pd.DataFrame()

    matches = []
    for _, row in lincs_compounds.iterrows():
        comp_name = str(row['compound']).lower().strip()

        # Direct name match
        foodb_match = foodb_compound[foodb_compound['name_lower'] == comp_name]

        # Try partial match (first word >= 4 chars)
        if foodb_match.empty:
            first_word = comp_name.split()[0] if comp_name else ""
            if len(first_word) >= 4:
                foodb_match = foodb_compound[
                    foodb_compound['name_lower'].str.contains(first_word, na=False)
                ].head(5)

        if not foodb_match.empty:
            for _, fmatch in foodb_match.iterrows():
                matches.append({
                    'lincs_compound': row['compound'],
                    'lincs_rank': row.get('rank', None),
                    'foodb_id': fmatch['id'],
                    'foodb_public_id': fmatch.get('public_id', ''),
                    'foodb_name': fmatch.get('name', ''),
                    'foodb_smiles': fmatch.get('moldb_smiles', ''),
                    'foodb_inchikey': fmatch.get('moldb_inchikey', ''),
                    'match_type': 'exact' if fmatch['name_lower'] == comp_name else 'partial'
                })

    return pd.DataFrame(matches)


def load_lincs_hits(stratum):
    """Load annotated LINCS reverser hits."""
    fpath = os.path.join(RESULTS_DIR, f"lincs_{stratum}", "tier1_ranking_annotated.csv")
    if not os.path.exists(fpath):
        fpath = os.path.join(RESULTS_DIR, f"lincs_{stratum}", "tier1_ranking.csv")
    if not os.path.exists(fpath):
        return pd.DataFrame()
    df = pd.read_csv(fpath)

    def is_valid(name):
        if not isinstance(name, str) or len(name) < 3:
            return False
        if re.match(r'^[\d,\s\.\-]+$', name):
            return False
        return True

    df = df[df['compound'].apply(is_valid)]
    return df


def main():
    print("=" * 70)
    print("  Nutrient-Gene Mapping: LINCS hits -> Food compounds -> Vegan sources")
    print("=" * 70)

    # Load FooDB
    print("\n--- Loading FooDB data ---")
    foodb_compound = load_foodb_compounds()
    foodb_food = load_foodb_foods()

    # Load consensus genes
    consensus_genes = {}
    for stratum in ['liver_immune', 'brain']:
        for direction in ['up', 'down']:
            fpath = os.path.join(PROC_DIR, f"signature_{stratum}_{direction}_genes.txt")
            if os.path.exists(fpath):
                with open(fpath) as f:
                    genes = {g.strip().upper() for g in f.read().strip().split('\n') if g.strip()}
                consensus_genes[f"{stratum}_{direction}"] = genes

    all_results = []

    for stratum in ['liver_immune', 'brain']:
        print(f"\n{'='*70}")
        print(f"  Processing: {stratum}")
        print(f"{'='*70}")

        lincs_df = load_lincs_hits(stratum)
        if lincs_df.empty:
            print(f"  No LINCS results for {stratum}")
            continue
        print(f"  LINCS hits: {len(lincs_df)} compounds")

        # Top 50 for nutrient mapping
        top_hits = lincs_df.head(50).copy()
        top_hits['rank'] = range(1, len(top_hits) + 1)

        # Match to FooDB
        print(f"  Matching compounds to FooDB...")
        compound_matches = match_compounds_to_foodb(top_hits, foodb_compound)
        n_matched = compound_matches['lincs_compound'].nunique() if not compound_matches.empty else 0
        print(f"  Matched {n_matched} compounds to FooDB")

        # Load content for matched compounds
        vegan_food_map = pd.DataFrame()
        if not compound_matches.empty:
            compound_ids = compound_matches['foodb_id'].unique()
            content = load_foodb_content(compound_ids)
            if content is not None and foodb_food is not None:
                # Merge with food names
                content = content.merge(
                    foodb_food[['id', 'name', 'food_group', 'food_subgroup']].add_prefix('food_'),
                    left_on='food_id', right_on='food_id',
                    how='left', suffixes=('', '_food')
                )
                # Also use orig_food_common_name if available
                content['food_name_full'] = content['food_name'].fillna(
                    content.get('orig_food_common_name', '')
                )
                content['is_vegan'] = content['food_name_full'].apply(is_vegan_food)
                vegan_food_map = content[content['is_vegan']].copy()
                print(f"  Vegan food-compound associations: {len(vegan_food_map)}")
                if not vegan_food_map.empty:
                    print(f"  Unique vegan foods: {vegan_food_map['food_name_full'].nunique()}")

        # Build nutrient-gene map
        target_col = [c for c in lincs_df.columns if 'drh_target' in c.lower()]
        moa_col = [c for c in lincs_df.columns if 'drh_moa' in c.lower()]
        phase_col = [c for c in lincs_df.columns if 'drh_clinical' in c.lower() or 'drh_phase' in c.lower()]

        stratum_genes = set()
        for key, genes in consensus_genes.items():
            if key.startswith(stratum):
                stratum_genes.update(genes)

        results = []
        for _, row in top_hits.iterrows():
            comp = row['compound']
            targets = row[target_col[0]] if target_col else None
            moa = row[moa_col[0]] if moa_col else None
            phase = row[phase_col[0]] if phase_col else None

            foodb_match = compound_matches[compound_matches['lincs_compound'] == comp] if not compound_matches.empty else pd.DataFrame()

            food_sources = []
            if not vegan_food_map.empty and not foodb_match.empty:
                foodb_ids = foodb_match['foodb_id'].tolist()
                foods = vegan_food_map[vegan_food_map['source_id'].isin(foodb_ids)]
                food_sources = foods['food_name_full'].dropna().unique().tolist()

            is_nutrient = comp.lower() in KNOWN_NUTRIENTS

            target_list = []
            if isinstance(targets, str) and targets != 'None':
                target_list = [t.strip() for t in re.split(r'[|;,]', targets)
                               if t.strip() and t.strip() != 'None']

            gene_overlap = [g for g in target_list if g in stratum_genes] if stratum_genes else []

            results.append({
                'stratum': stratum,
                'compound': comp,
                'lincs_rank': row.get('rank', None),
                'best_z_sum': row.get('best_z_sum', None),
                'n_reversing_sigs': row.get('n_reversing_sigs', None),
                'moa': moa if isinstance(moa, str) and moa != 'None' else None,
                'clinical_phase': phase if isinstance(phase, str) and phase != 'None' else None,
                'drug_targets': '; '.join(target_list) if target_list else None,
                'n_targets': len(target_list),
                'signature_gene_overlap': '; '.join(gene_overlap) if gene_overlap else None,
                'n_gene_overlap': len(gene_overlap),
                'foodb_matched': not foodb_match.empty,
                'foodb_name': foodb_match['foodb_name'].iloc[0] if not foodb_match.empty else None,
                'is_known_nutrient': is_nutrient,
                'vegan_food_sources': '; '.join(food_sources[:15]) if food_sources else None,
                'n_vegan_foods': len(food_sources),
            })

        nutrient_map = pd.DataFrame(results)
        if not nutrient_map.empty:
            nutrient_map.to_csv(
                os.path.join(RESULTS_DIR, f"nutrient_gene_map_{stratum}.csv"),
                index=False
            )

            n_with_food = nutrient_map['foodb_matched'].sum()
            n_nutrients = nutrient_map['is_known_nutrient'].sum()
            n_with_vegan = (nutrient_map['n_vegan_foods'] > 0).sum()
            n_gene_overlap = (nutrient_map['n_gene_overlap'] > 0).sum()

            print(f"\n  Summary for {stratum}:")
            print(f"    Compounds analyzed: {len(nutrient_map)}")
            print(f"    Matched to FooDB: {n_with_food}")
            print(f"    Known nutrients: {n_nutrients}")
            print(f"    With vegan food sources: {n_with_vegan}")
            print(f"    With signature gene overlap: {n_gene_overlap}")

            # Print nutrient hits
            nutrient_hits = nutrient_map[
                nutrient_map['is_known_nutrient'] | (nutrient_map['n_vegan_foods'] > 0)
            ]
            if not nutrient_hits.empty:
                print(f"\n  Top nutrient/food-compound hits:")
                display_cols = ['compound', 'lincs_rank', 'moa', 'is_known_nutrient',
                                'n_vegan_foods', 'vegan_food_sources']
                available = [c for c in display_cols if c in nutrient_hits.columns]
                for _, r in nutrient_hits.head(15).iterrows():
                    print(f"    #{int(r['lincs_rank'])} {r['compound']}"
                          f" | nutrient={r['is_known_nutrient']}"
                          f" | vegan_foods={r['n_vegan_foods']}"
                          f" | moa={r.get('moa', 'N/A')}")

            all_results.append(nutrient_map)

    # Combined
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.to_csv(os.path.join(RESULTS_DIR, "nutrient_gene_map_combined.csv"), index=False)
        print(f"\n{'='*70}")
        print(f"  Combined nutrient-gene map: {len(combined)} entries")
        print(f"{'='*70}")

    print("\nNext step: Run 06_recipe_optimization.py for vegan recipe generation.")


if __name__ == "__main__":
    main()
