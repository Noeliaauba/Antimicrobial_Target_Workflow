#!/usr/bin/env python


from pathlib import Path
import argparse

def extract_tm_ids(predicted_file: Path, min_tm: int = 2):
    tm_ids = set()

    with predicted_file.open() as infile:
        next(infile)  # skip header

        for line in infile:
            line = line.strip()
            if not line:
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            protein_id = parts[0]

            try:
                num_tm = int(parts[1])
                sp_flag = int(parts[2])  # 🔥 ahora es int
            except ValueError:
                continue

            # 🔥 TU CONDICIÓN EXACTA
            if num_tm >=1 and sp_flag == 0:
                tm_ids.add(protein_id)

    return tm_ids

# Write Fasta enries whose ID is in tm_ids
def filter_by_id(faa_file: Path, output_file: Path, valid_ids):
    written = 0
    total = 0

    with faa_file.open() as infile, output_file.open("w") as outfile:
        write = False

        for line in infile:
            if line.startswith(">"):
                total += 1
                protein_id = line.split()[0].lstrip(">")
                write = protein_id in valid_ids

                if write:
                    written += 1
                    

            if write:
                outfile.write(line)

    print(f"{faa_file.name}:")
    print(f"   total: {total}")
    print(f"   TM proteins: {written}")
    print(f"   percentage: {written/total:.2%}\n")

    return total, written

def process_paths(plasmid_dir: Path, topology_dir: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    plasmid_files = list(plasmid_dir.glob("*.faa"))
    total_all = 0
    tm_all = 0

    for p_file in plasmid_files:
        sample_name = p_file.stem
        topology_file = topology_dir / f"{sample_name}.txt"

        if not topology_file.exists():
            print(f"[WARNING] No topology file for {sample_name}")
            continue

        tm_ids = extract_tm_ids(topology_file, min_tm=1)

        output_file = output_dir / p_file.name
        total, written = filter_by_id(p_file, output_file, tm_ids)

        total_all += total
        tm_all += written

    # GLOBAL SUMMARY
    percentage_all = (tm_all / total_all) if total_all > 0 else 0

    print("\nGLOBAL SUMMARY")
    print(f"Total proteins: {total_all}")
    print(f"TM proteins: {tm_all}")
    print(f"Percentage: {percentage_all:.2%}")



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plasmid_dir", type=Path, required=True)
    parser.add_argument("--topology_dir", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    process_paths(args.plasmid_dir, args.topology_dir, args.output_dir)


if __name__ == "__main__":
    main()











