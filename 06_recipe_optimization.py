#!/usr/bin/env python3
"""
06_recipe_optimization.py — Vegan recipe optimization for spaceflight
countermeasure nutrients.

Since RecipeNLG requires manual download (not programmatically accessible),
this script builds a targeted vegan recipe database from:
1. FooDB food-compound associations (from Script 05)
2. A curated vegan recipe template library targeting key nutrients
3. Phase-based scoring (pre-flight, in-flight, post-flight)

The script generates optimized vegan meal plans that maximize coverage of
spaceflight countermeasure nutrients (quercetin, sulforaphane, curcumin,
ascorbic acid, luteolin, etc.) identified in Script 05.
"""
import json, os, re
import pandas as pd
import numpy as np
from collections import defaultdict

PROJECT_ROOT = "/workspace/astronaut-opposite-forcing"
PROC_DIR     = os.path.join(PROJECT_ROOT, "data", "processed")
RESULTS_DIR  = os.path.join(PROJECT_ROOT, "results", "tables")
FIGURES_DIR  = os.path.join(PROJECT_ROOT, "results", "figures")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Key spaceflight countermeasure nutrients from LINCS + FooDB ────────────
COUNTERMEASURE_NUTRIENTS = {
    # Nutrient: (description, primary vegan sources, target pathways)
    "quercetin": {
        "description": "Flavonoid antioxidant; senolytic; AMPK activator",
        "key_sources": ["onion", "apple", "broccoli", "kale", "capers", "berry"],
        "pathways": ["Nrf2/ARE", "AMPK", "senescence clearance"],
        "lincs_stratum": "liver_immune",
    },
    "sulforaphane": {
        "description": "Isothiocyanate; Nrf2 activator; anti-inflammatory",
        "key_sources": ["broccoli", "brussels sprouts", "cabbage", "kale"],
        "pathways": ["Nrf2/ARE", "Phase II detox", "NF-kB inhibition"],
        "lincs_stratum": "brain",
    },
    "curcumin": {
        "description": "Polyphenol; anti-inflammatory; NF-kB inhibitor",
        "key_sources": ["turmeric", "curry"],
        "pathways": ["NF-kB", "Nrf2/ARE", "COX-2 inhibition"],
        "lincs_stratum": "brain",
    },
    "ascorbic acid": {
        "description": "Vitamin C; antioxidant; collagen synthesis",
        "key_sources": ["orange", "lemon", "pepper", "broccoli", "kiwi", "strawberry"],
        "pathways": ["antioxidant", "collagen synthesis", "immune support"],
        "lincs_stratum": "brain",
    },
    "luteolin": {
        "description": "Flavone; antioxidant; anti-inflammatory",
        "key_sources": ["celery", "parsley", "pepper", "oregano", "thyme"],
        "pathways": ["Nrf2/ARE", "NF-kB inhibition", "MAO inhibition"],
        "lincs_stratum": "brain",
    },
    "resveratrol": {
        "description": "Stilbene; SIRT1 activator; mitochondrial biogenesis",
        "key_sources": ["grape", "peanut", "berry"],
        "pathways": ["SIRT1", "AMPK", "mitochondrial biogenesis"],
        "lincs_stratum": "liver_immune",
    },
    "epigallocatechin gallate": {
        "description": "Catechin; antioxidant; AMPK activator",
        "key_sources": ["green tea"],
        "pathways": ["AMPK", "Nrf2/ARE", "mTOR inhibition"],
        "lincs_stratum": "liver_immune",
    },
    "gingerol": {
        "description": "Phenol; anti-inflammatory; anti-nausea",
        "key_sources": ["ginger"],
        "pathways": ["NF-kB", "COX-2 inhibition", "5-HT3 antagonism"],
        "lincs_stratum": "liver_immune",
    },
    "allicin": {
        "description": "Organosulfur; antimicrobial; cardiovascular",
        "key_sources": ["garlic", "onion"],
        "pathways": ["cardiovascular", "antimicrobial", "antioxidant"],
        "lincs_stratum": "liver_immune",
    },
    "lutein": {
        "description": "Carotenoid; eye health; antioxidant",
        "key_sources": ["kale", "spinach", "broccoli"],
        "pathways": ["antioxidant", "macular protection"],
        "lincs_stratum": "brain",
    },
}

# ── Curated vegan recipe templates ─────────────────────────────────────────
# Each recipe targets specific countermeasure nutrients
VEGAN_RECIPES = [
    {
        "name": "Broccoli Sulforaphane Power Bowl",
        "meal_type": "lunch",
        "phase": "in_flight",
        "ingredients": [
            "broccoli (raw or lightly steamed)", "quinoa", "avocado",
            "sunflower seeds", "lemon juice", "garlic", "olive oil",
            "turmeric", "black pepper"
        ],
        "target_nutrients": ["sulforaphane", "curcumin", "ascorbic acid", "quercetin", "lutein"],
        "preparation": "Lightly steam broccoli (maximizes sulforaphane). Cook quinoa. "
                       "Toss with avocado, seeds, lemon-garlic-turmeric dressing.",
        "nutrient_score": 0,
        "countermeasure_tags": ["Nrf2 activation", "antioxidant", "anti-inflammatory"],
    },
    {
        "name": "Golden Turmeric Lentil Soup",
        "meal_type": "dinner",
        "phase": "in_flight",
        "ingredients": [
            "red lentils", "turmeric", "ginger", "garlic", "onion",
            "coconut milk", "spinach", "lemon juice", "cumin", "coriander"
        ],
        "target_nutrients": ["curcumin", "gingerol", "allicin", "ascorbic acid", "quercetin"],
        "preparation": "Sauté onion, garlic, ginger. Add turmeric, cumin, coriander. "
                       "Add lentils and water, cook 20 min. Stir in spinach and coconut milk.",
        "nutrient_score": 0,
        "countermeasure_tags": ["anti-inflammatory", "NF-kB inhibition", "protein-rich"],
    },
    {
        "name": "Green Tea Berry Antioxidant Smoothie",
        "meal_type": "breakfast",
        "phase": "in_flight",
        "ingredients": [
            "green tea (brewed, chilled)", "blueberry", "strawberry",
            "banana", "flax seed", "spinach", "lemon juice"
        ],
        "target_nutrients": ["epigallocatechin gallate", "ascorbic acid", "quercetin", "lutein"],
        "preparation": "Brew green tea, chill. Blend with berries, banana, spinach, flax, lemon.",
        "nutrient_score": 0,
        "countermeasure_tags": ["AMPK activation", "antioxidant", "mitochondrial support"],
    },
    {
        "name": "Kale Caesar with Capers and Nutritional Yeast",
        "meal_type": "lunch",
        "phase": "pre_flight",
        "ingredients": [
            "kale (massaged)", "capers", "nutritional yeast", "garlic",
            "lemon juice", "olive oil", "hemp seed", "mustard", "black pepper"
        ],
        "target_nutrients": ["quercetin", "luteolin", "ascorbic acid", "lutein", "allicin"],
        "preparation": "Massage kale with olive oil and lemon. Add capers, nutritional yeast, "
                       "hemp seeds. Toss with garlic-mustard dressing.",
        "nutrient_score": 0,
        "countermeasure_tags": ["senolytic", "Nrf2 activation", "B12-fortified"],
    },
    {
        "name": "Spiced Chickpea and Spinach Curry",
        "meal_type": "dinner",
        "phase": "in_flight",
        "ingredients": [
            "chickpeas", "spinach", "tomato", "onion", "garlic", "ginger",
            "turmeric", "cumin", "coriander", "cinnamon", "coconut milk"
        ],
        "target_nutrients": ["curcumin", "gingerol", "allicin", "quercetin", "lutein", "ascorbic acid"],
        "preparation": "Sauté onion, garlic, ginger with spices. Add tomatoes, chickpeas, "
                       "coconut milk. Simmer 15 min. Add spinach at end.",
        "nutrient_score": 0,
        "countermeasure_tags": ["anti-inflammatory", "protein-rich", "iron-rich"],
    },
    {
        "name": "Cruciferous Crunch Salad with Tahini",
        "meal_type": "lunch",
        "phase": "post_flight",
        "ingredients": [
            "broccoli (raw, grated)", "brussels sprouts (shaved)", "cabbage",
            "carrot", "tahini", "lemon juice", "garlic", "apple cider vinegar",
            "pomegranate seeds", "walnut"
        ],
        "target_nutrients": ["sulforaphane", "quercetin", "ascorbic acid", "allicin"],
        "preparation": "Grate broccoli and brussels sprouts. Shave cabbage. "
                       "Toss with tahini-lemon-garlic dressing. Top with pomegranate and walnuts.",
        "nutrient_score": 0,
        "countermeasure_tags": ["Nrf2 activation", "senolytic", "antioxidant"],
    },
    {
        "name": "Ginger-Green Tea Recovery Broth",
        "meal_type": "snack",
        "phase": "post_flight",
        "ingredients": [
            "ginger (fresh)", "green tea", "lemon", "turmeric",
            "black pepper", "honey substitute (agave)", "seaweed (nori)"
        ],
        "target_nutrients": ["gingerol", "epigallocatechin gallate", "curcumin", "ascorbic acid"],
        "preparation": "Brew ginger and green tea. Add turmeric, black pepper, lemon, agave. "
                       "Float nori strips.",
        "nutrient_score": 0,
        "countermeasure_tags": ["anti-nausea", "mitochondrial support", "recovery"],
    },
    {
        "name": "Overnight Oats with Berries and Flax",
        "meal_type": "breakfast",
        "phase": "pre_flight",
        "ingredients": [
            "rolled oats", "flax seed", "chia seed", "blueberry",
            "strawberry", "walnut", "cinnamon", "vanilla", "soy milk"
        ],
        "target_nutrients": ["quercetin", "ascorbic acid", "epigallocatechin gallate"],
        "preparation": "Combine oats, seeds, soy milk, cinnamon, vanilla. Refrigerate overnight. "
                       "Top with berries and walnuts in the morning.",
        "nutrient_score": 0,
        "countermeasure_tags": ["omega-3", "fiber-rich", "antioxidant"],
    },
    {
        "name": "Mediterranean Quinoa with Herbs",
        "meal_type": "lunch",
        "phase": "in_flight",
        "ingredients": [
            "quinoa", "parsley", "oregano", "thyme", "tomato", "cucumber",
            "olive oil", "lemon juice", "olive", "red onion", "mint"
        ],
        "target_nutrients": ["luteolin", "quercetin", "ascorbic acid", "allicin"],
        "preparation": "Cook quinoa. Chop herbs and vegetables. Toss with olive oil-lemon dressing.",
        "nutrient_score": 0,
        "countermeasure_tags": ["anti-inflammatory", "Mediterranean diet", "antioxidant"],
    },
    {
        "name": "Roasted Root Vegetables with Rosemary",
        "meal_type": "dinner",
        "phase": "pre_flight",
        "ingredients": [
            "sweet potato", "carrot", "beet", "onion", "garlic",
            "rosemary", "thyme", "olive oil", "black pepper", "balsamic vinegar"
        ],
        "target_nutrients": ["quercetin", "luteolin", "allicin", "ascorbic acid"],
        "preparation": "Cube vegetables. Toss with olive oil, herbs, garlic. Roast at 200°C "
                       "for 35 min. Drizzle with balsamic.",
        "nutrient_score": 0,
        "countermeasure_tags": ["beta-carotene", "antioxidant", "fiber-rich"],
    },
    {
        "name": "Spicy Peanut Soba Noodles",
        "meal_type": "dinner",
        "phase": "in_flight",
        "ingredients": [
            "buckwheat soba noodles", "peanut butter", "ginger", "garlic",
            "soy sauce", "rice vinegar", "chili", "cabbage", "carrot", "scallion"
        ],
        "target_nutrients": ["resveratrol", "gingerol", "allicin", "quercetin"],
        "preparation": "Cook soba noodles. Whisk peanut sauce with ginger, garlic, soy, vinegar. "
                       "Toss with shredded cabbage and carrot.",
        "nutrient_score": 0,
        "countermeasure_tags": ["protein-rich", "anti-inflammatory", "resveratrol"],
    },
    {
        "name": "Cocoa-Chia Pudding with Raspberries",
        "meal_type": "breakfast",
        "phase": "post_flight",
        "ingredients": [
            "chia seed", "cocoa", "raspberry", "almond milk",
            "vanilla", "maple syrup", "walnut"
        ],
        "target_nutrients": ["quercetin", "ascorbic acid"],
        "preparation": "Mix chia, cocoa, almond milk, vanilla, maple syrup. Refrigerate 4+ hours. "
                       "Top with raspberries and walnuts.",
        "nutrient_score": 0,
        "countermeasure_tags": ["omega-3", "antioxidant", "flavonoid-rich"],
    },
]


def score_recipe(recipe, nutrient_map_df):
    """Score a recipe based on nutrient coverage and countermeasure relevance."""
    score = 0
    nutrients_hit = []

    for nutrient in recipe["target_nutrients"]:
        if nutrient in COUNTERMEASURE_NUTRIENTS:
            # Base score for each countermeasure nutrient
            score += 10

            # Bonus if nutrient was identified in LINCS screening
            nutrient_info = COUNTERMEASURE_NUTRIENTS[nutrient]
            stratum = nutrient_info["lincs_stratum"]

            # Check if this nutrient appeared in LINCS results
            if not nutrient_map_df.empty:
                nutrient_in_lincs = nutrient_map_df[
                    nutrient_map_df['compound'].str.lower().str.contains(
                        nutrient.split()[0], na=False
                    )
                ]
                if not nutrient_in_lincs.empty:
                    score += 15  # Bonus for LINCS-validated nutrient

            # Bonus for number of vegan food sources
            n_sources = len(nutrient_info["key_sources"])
            score += min(n_sources, 5)

            nutrients_hit.append(nutrient)

    # Bonus for number of countermeasure tags
    score += len(recipe["countermeasure_tags"]) * 3

    # Phase relevance: in_flight recipes get bonus (most critical phase)
    if recipe["phase"] == "in_flight":
        score += 5
    elif recipe["phase"] == "pre_flight":
        score += 3

    return score, nutrients_hit


def build_meal_plans(scored_recipes):
    """Build phase-based meal plans from scored recipes."""
    phases = ["pre_flight", "in_flight", "post_flight"]
    meal_types = ["breakfast", "lunch", "dinner", "snack"]

    plans = {}
    for phase in phases:
        phase_recipes = [r for r in scored_recipes if r["phase"] == phase]
        phase_recipes.sort(key=lambda x: x["nutrient_score"], reverse=True)

        plan = {
            "phase": phase,
            "recipes": [],
            "total_nutrients": set(),
            "avg_score": 0,
        }

        # Select top recipes for each meal type
        for mt in meal_types:
            mt_recipes = [r for r in phase_recipes if r["meal_type"] == mt]
            if mt_recipes:
                best = mt_recipes[0]
                plan["recipes"].append({
                    "meal_type": mt,
                    "name": best["name"],
                    "ingredients": best["ingredients"],
                    "preparation": best["preparation"],
                    "target_nutrients": best["target_nutrients"],
                    "nutrient_score": best["nutrient_score"],
                    "countermeasure_tags": best["countermeasure_tags"],
                })
                plan["total_nutrients"].update(best["target_nutrients"])

        plan["total_nutrients"] = sorted(plan["total_nutrients"])
        plan["avg_score"] = np.mean([r["nutrient_score"] for r in phase_recipes]) if phase_recipes else 0
        plans[phase] = plan

    return plans


def main():
    print("=" * 70)
    print("  Vegan Recipe Optimization for Spaceflight Countermeasures")
    print("=" * 70)

    # Load nutrient-gene map
    nutrient_map_path = os.path.join(RESULTS_DIR, "nutrient_gene_map_combined.csv")
    nutrient_map_df = pd.DataFrame()
    if os.path.exists(nutrient_map_path):
        nutrient_map_df = pd.read_csv(nutrient_map_path)
        print(f"  Loaded nutrient-gene map: {len(nutrient_map_df)} entries")
    else:
        print(f"  Nutrient-gene map not found — using nutrient knowledge base only")

    # Score recipes
    print(f"\n--- Scoring {len(VEGAN_RECIPES)} vegan recipes ---")
    scored_recipes = []
    for recipe in VEGAN_RECIPES:
        score, nutrients_hit = score_recipe(recipe, nutrient_map_df)
        recipe["nutrient_score"] = score
        recipe["nutrients_hit"] = nutrients_hit
        scored_recipes.append(recipe)
        print(f"  {recipe['name']}: score={score}, nutrients={nutrients_hit}")

    # Sort by score
    scored_recipes.sort(key=lambda x: x["nutrient_score"], reverse=True)

    # Build meal plans
    print(f"\n--- Building phase-based meal plans ---")
    meal_plans = build_meal_plans(scored_recipes)

    for phase, plan in meal_plans.items():
        print(f"\n  {phase.upper()} ({len(plan['recipes'])} meals, avg score={plan['avg_score']:.1f}):")
        print(f"    Nutrients covered: {', '.join(plan['total_nutrients'])}")
        for meal in plan["recipes"]:
            print(f"    [{meal['meal_type']}] {meal['name']} (score={meal['nutrient_score']})")

    # Save outputs
    # 1. Scored recipes
    recipes_df = pd.DataFrame([
        {
            "name": r["name"],
            "meal_type": r["meal_type"],
            "phase": r["phase"],
            "ingredients": "; ".join(r["ingredients"]),
            "target_nutrients": "; ".join(r["target_nutrients"]),
            "nutrient_score": r["nutrient_score"],
            "countermeasure_tags": "; ".join(r["countermeasure_tags"]),
            "preparation": r["preparation"],
        }
        for r in scored_recipes
    ])
    recipes_df.to_csv(os.path.join(RESULTS_DIR, "vegan_recipes_scored.csv"), index=False)

    # 2. Meal plans
    plan_records = []
    for phase, plan in meal_plans.items():
        for meal in plan["recipes"]:
            plan_records.append({
                "phase": phase,
                "meal_type": meal["meal_type"],
                "recipe_name": meal["name"],
                "ingredients": "; ".join(meal["ingredients"]),
                "preparation": meal["preparation"],
                "target_nutrients": "; ".join(meal["target_nutrients"]),
                "nutrient_score": meal["nutrient_score"],
                "countermeasure_tags": "; ".join(meal["countermeasure_tags"]),
            })
    plans_df = pd.DataFrame(plan_records)
    plans_df.to_csv(os.path.join(RESULTS_DIR, "vegan_meal_plans.csv"), index=False)

    # 3. Nutrient summary
    nutrient_summary = []
    for nutrient, info in COUNTERMEASURE_NUTRIENTS.items():
        recipes_with = [r["name"] for r in scored_recipes if nutrient in r["target_nutrients"]]
        nutrient_summary.append({
            "nutrient": nutrient,
            "description": info["description"],
            "key_vegan_sources": "; ".join(info["key_sources"]),
            "target_pathways": "; ".join(info["pathways"]),
            "lincs_stratum": info["lincs_stratum"],
            "n_recipes_targeting": len(recipes_with),
            "recipes": "; ".join(recipes_with),
        })
    nutrient_summary_df = pd.DataFrame(nutrient_summary)
    nutrient_summary_df.to_csv(os.path.join(RESULTS_DIR, "countermeasure_nutrient_summary.csv"), index=False)

    # 4. Full JSON output
    output = {
        "recipes": [
            {
                "name": r["name"],
                "meal_type": r["meal_type"],
                "phase": r["phase"],
                "ingredients": r["ingredients"],
                "preparation": r["preparation"],
                "target_nutrients": r["target_nutrients"],
                "nutrient_score": r["nutrient_score"],
                "countermeasure_tags": r["countermeasure_tags"],
            }
            for r in scored_recipes
        ],
        "meal_plans": {
            phase: {
                "phase": plan["phase"],
                "total_nutrients": plan["total_nutrients"],
                "avg_score": plan["avg_score"],
                "recipes": plan["recipes"],
            }
            for phase, plan in meal_plans.items()
        },
        "countermeasure_nutrients": {
            k: {
                "description": v["description"],
                "key_sources": v["key_sources"],
                "pathways": v["pathways"],
                "lincs_stratum": v["lincs_stratum"],
            }
            for k, v in COUNTERMEASURE_NUTRIENTS.items()
        },
    }
    with open(os.path.join(RESULTS_DIR, "vegan_recipe_optimization.json"), "w") as f:
        json.dump(output, f, indent=2)

    print(f"\n{'='*70}")
    print(f"  Recipe optimization complete!")
    print(f"  {len(scored_recipes)} recipes scored, {len(meal_plans)} phase plans built")
    print(f"  {len(COUNTERMEASURE_NUTRIENTS)} countermeasure nutrients targeted")
    print(f"  Results saved to: {RESULTS_DIR}/")
    print(f"{'='*70}")

    # Print top 5 recipes
    print(f"\n  Top 5 recipes by nutrient score:")
    for r in scored_recipes[:5]:
        print(f"    [{r['nutrient_score']}] {r['name']} ({r['phase']})")
        print(f"       Nutrients: {', '.join(r['target_nutrients'])}")
        print(f"       Tags: {', '.join(r['countermeasure_tags'])}")

    print("\nNext step: Run 07_report_generation.py for report + Zenodo packaging.")


if __name__ == "__main__":
    main()
