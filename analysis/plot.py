#Standard libraries
import argparse
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg") # Generate plots without requiring a graphical display
import matplotlib.pyplot as plt
import numpy as np      

# Load JSONL statistics file
def load_jsonl(path: Path):
    data = []
    with open(path) as f:
        for line in f:
            # Extract each line as a JSON object into the data list
            record = json.loads(line)
            data.append(record)
    return data
       

# Generation of barplot
def plot_barplot(data, title, output_dir, base_name):
    output_path = (output_dir / f"{base_name}_barplot.png")
    # Extract these values from statistics file
    names = [d["genome"] for d in data]
    total = [d["total"] for d in data]
    kept = [d["kept"] for d in data]
    # Calculate the global counts
    total_sum = sum(total)
    kept_sum = sum(kept)
    # Configuration of the bars
    x = np.arange(len(names))
    width = 0.25
    offset = 0.01
    plt.figure(figsize=(12, 6))
    # Creation of the "total proteins" bars and the "filtered proteins" bars 
    total_bars = plt.bar(x - width/2, total, width, label=f"Total ({total_sum:,})")
    kept_bars = plt.bar(x + width/2, kept, width, label=f"Kept ({kept_sum:,})")
    # Configuration of the labels on top of the bars with the counts of total and kept proteins
    for bar in total_bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, height + 5, f"{int(height)}",ha="center", va="bottom",fontsize=10)
    for bar in kept_bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2+offset, height + 5, f"{int(height)}", ha="left", va="bottom",fontsize=10)
    # Configuration of the barplot 
    plt.xticks(x, names, rotation=90)
    plt.ylabel("Proteins")
    plt.xlabel("Genomes")
    plt.title(f"{title} Effect")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


# Generation of heatmap
def plot_heatmap(data,title, output_dir, base_name):
    output_path = (output_dir / f"{base_name}_heatmap.png")
    # Extract these values from statistics file
    names = [d["genome"] for d in data]
    percentages = np.array([[d["percentage"]] for d in data])
    # Creates heatmap
    plt.figure(figsize=(5, 6))
    im = plt.imshow(percentages, aspect="auto" )
    # Configuration of axes and colorbar
    plt.yticks(range(len(names)),names, fontsize=8)
    plt.xticks([0],["Retention %"])
    plt.colorbar(im,label="Percentage")
    plt.title(f"{title} Retention")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()


# Parse the input arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Generate filtering plots" )
    parser.add_argument("--stats_dir", type=Path, required=True)
    parser.add_argument("--plots_dir", type=Path, required=True)
    return parser.parse_args()

# Retrieve all JSONL statistics files excluding cluster statistics 
def get_stats_files(stats_dir: Path):
    stats_files = []
    for file in sorted(stats_dir.glob("*.jsonl")):
        if file.name == "clusters.jsonl":
            continue
        stats_files.append(file)

    return stats_files

def main():
    # Read the parameters
    args = parse_args()
    stats_dir = args.stats_dir
    plots_dir = args.plots_dir
    plots_dir.mkdir(parents=True, exist_ok=True)
    # Retrieve statistics files
    stats_files = get_stats_files(stats_dir)
    if not stats_files:
        print("[ERROR] No JSONL files found.")
        return
    # Generate both plots for each filtering step
    for stats_file in stats_files:
        print(f"Processing {stats_file.name}")
        data = load_jsonl(stats_file)
        base_name = stats_file.stem
        plot_barplot(data,base_name,plots_dir,base_name)
        plot_heatmap(data,base_name, plots_dir, base_name)

if __name__ == "__main__":
    main()