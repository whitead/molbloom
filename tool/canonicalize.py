import os
import multiprocessing as mp
from pathlib import Path
import time
from tqdm import tqdm
from molbloom import canon


def process_smiles_file(file_path):
    input_path = Path(file_path)
    output_path = input_path.with_name(f"{input_path.stem}_canonical.smi")

    with open(input_path, "r") as f:
        smiles_data = f.readlines()

    canonical_smiles = [
        canon(s.split()[0].strip()) for s in smiles_data[1:] if s.strip()
    ]

    with open(output_path, "w") as f:
        f.writelines([c + "\n" for c in canonical_smiles if c])

    return (str(input_path), str(output_path), True)


def apply_canonical_transformation(smiles_data):
    # Replace with your actual transformation logic
    transformed_data = smiles_data
    return transformed_data


def process_smiles_files(file_list, num_processes=None):
    if num_processes is None:
        num_processes = mp.cpu_count()

    pool = mp.Pool(processes=num_processes)
    start_time = time.time()
    results = []

    with tqdm(total=len(file_list), desc="Processing SMILES files") as pbar:
        for result in pool.imap_unordered(process_smiles_file, file_list):
            results.append(result)
            pbar.update(1)

    pool.close()
    pool.join()

    end_time = time.time()
    successful = sum(1 for _, _, success in results if success)

    print(
        f"\nComplete: {successful}/{len(file_list)} files processed in {end_time - start_time:.2f}s"
    )

    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Process SMILES files to canonical form."
    )
    parser.add_argument("files", nargs="+", help="SMILES files to process")
    parser.add_argument(
        "--processes", type=int, default=None, help="Number of processes"
    )
    args = parser.parse_args()

    valid_files = []
    for pattern in args.files:
        if os.path.isfile(pattern) and pattern.endswith(".smi"):
            valid_files.append(pattern)
        elif "*" in pattern:
            from glob import glob

            matches = glob(pattern)
            valid_files.extend(
                [f for f in matches if os.path.isfile(f) and f.endswith(".smi")]
            )

    # remove those with canonical in their name
    valid_files = [f for f in valid_files if "_canonical" not in f]
    if not valid_files:
        print("No valid .smi files to process after filtering!")
        exit(1)

    print(f"Found {len(valid_files)} SMILES files to process")
    process_smiles_files(valid_files, args.processes)
