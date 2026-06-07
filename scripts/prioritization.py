# Standard libraries
import argparse
from pathlib import Path
import csv

#Local imports form utils.py to reuse declared functions
from utils import ensure_dir

GO_SCORES = {
    "transport": {"terms": {"GO:0055085", "GO:0005215"},"score": 3},
    "ion_transport": {"terms": {"GO:0015075"},"score": 3},
    "membrane_localization": {"terms": {"GO:0016020", "GO:0005886", "GO:0071944"}, "score": 2 }
}

KO_SCORES = {
    "adaptation": {"terms": {"ko02024", "ko02020"}, "score": 3 },
    "nutrient_uptake": {"terms": {"ko02060"},"score": 3},
    "secretion": {"terms": {"ko03060", "ko03070"}, "score": 2}
}

MAP_SCORES = {
    "ABC_transport": {"terms": {"map02010"}, "score": 3},
    "adaptation": {"terms": {"map02024","map02020","map02060"},"score": 3},
    "secretion": {"terms": { "map03060","map03070"},"score": 2},
    "cell_wall": {"terms": {"map00550"}, "score": 2},
    "energy": {"terms": {"map00190"}, "score": 1},
}

# Separate terms 
def split_terms(value):
    if not value or value == "-":
        return []
    return [x.strip() for x in value.split(",")]


# Calculate score for matching categories
def score_categories(terms, categories):
    score = 0
    matched = []
    term_set = set(terms)
    for category, config in categories.items():
        if term_set.intersection(config["terms"]):
            score += config["score"]
            matched.append(category)
    return score, matched


def prioritize(input_file: Path, output_file: Path, go_threshold: int):
    results = []
    with open(input_file) as infile:
        lines = []
        for line in infile:
            if line.startswith("##"):
                continue
            if line.startswith("#query"): 
                line = line[1:] # Removes #, converting "#query" -> "query"
            lines.append(line)
        reader = csv.DictReader(lines, delimiter="\t")
        for row in reader:
            # GO terms
            go_terms = split_terms(row.get("GOs", ""))
            go_score, go_hits = score_categories(go_terms,GO_SCORES)
            # GO pre-filter
            if go_score < go_threshold:
                continue
            # KEGG terms (ko and map)
            kegg_terms = split_terms(row.get("KEGG_Pathway", ""))
            ko_terms = []
            for term in kegg_terms:
                term = term.strip()
                if term.startswith("ko"):
                    ko_terms.append(term)
            map_terms = []        
            for term in kegg_terms:
                term = term.strip()
                if term.startswith("map"):
                    map_terms.append(term)

            ko_score, ko_hits = score_categories(ko_terms, KO_SCORES)
            map_score, map_hits = score_categories(map_terms, MAP_SCORES)
            # Scores
            total_score = go_score + ko_score + map_score
            results.append({"protein_id": row["query"],"go_score": go_score,"ko_score": ko_score,"map_score": map_score,"total_score": total_score})
    results.sort(key=lambda x: x["total_score"],reverse=True)
    max_score = results[0]["total_score"]
    high_threshold = max_score * 0.7
    medium_threshold = max_score * 0.4
    print(f"[INFO] HIGH threshold: {high_threshold:.1f}")
    for row in results:
        if row["total_score"] >= high_threshold:
            row["priority"] = "HIGH"
        elif row["total_score"] >= medium_threshold:
            row["priority"] = "MEDIUM"
        else:
            row["priority"] = "LOW"
    ensure_dir(output_file.parent)
    with open(output_file, "w", newline="") as outfile:
        writer = csv.DictWriter(outfile, fieldnames=["protein_id","go_score","ko_score","map_score","total_score","priority"],delimiter="\t")
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    high_candidates = sum(1 for r in results if r["priority"] == "HIGH")


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser( description="Prioritization of antimicrobial target candidates")
    parser.add_argument("--input_file",type=Path,required=True)
    parser.add_argument("--output_file",type=Path,required=True)
    parser.add_argument("--go_threshold", type=int, default=5)
    return parser.parse_args()


def main():
    # Read the parameters
    args = parse_args()
    prioritize(args.input_file, args.output_file, args.go_threshold)

if __name__ == "__main__":
    main()

