"""Build a SureChEMBL bloom filter from a parquet file of compounds."""

import argparse
import multiprocessing as mp
import math

import pyarrow.parquet as pq
from tqdm import tqdm

from molbloom import canon
from molbloom.bloom import CustomFilter


FILTER_SIZE_BITS = 200 * 1024 * 1024 * 8  # 200 MB
BATCH_SIZE = 100_000


def find_smiles_column(schema):
    """Auto-detect the SMILES column from the parquet schema."""
    candidates = [f.name for f in schema if "smiles" in f.name.lower()]
    if len(candidates) == 1:
        return candidates[0]
    print("Schema:")
    for f in schema:
        print(f"  {f.name}: {f.type}")
    if len(candidates) == 0:
        raise SystemExit("Error: no column containing 'smiles' found.")
    raise SystemExit(f"Error: ambiguous SMILES columns: {candidates}")


def canon_batch(smiles_list):
    """Canonicalize a batch of SMILES strings."""
    return [canon(s) for s in smiles_list]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parquet", help="Path to compounds.parquet")
    parser.add_argument(
        "--output", default="sure-chembl.bloom", help="Output bloom file"
    )
    parser.add_argument(
        "--processes", type=int, default=None, help="Number of worker processes"
    )
    args = parser.parse_args()

    pf = pq.ParquetFile(args.parquet)
    schema = pf.schema_arrow
    smiles_col = find_smiles_column(schema)
    print(f"Using SMILES column: {smiles_col}")

    num_rows = pf.metadata.num_rows
    print(f"Total rows: {num_rows:,}")

    bf = CustomFilter(FILTER_SIZE_BITS, num_rows, "sure-chembl")

    nprocs = args.processes or max(1, mp.cpu_count() - 2)
    print(f"Using {nprocs} worker processes")

    added = 0
    failed = 0

    with mp.Pool(processes=nprocs) as pool:
        with tqdm(total=num_rows, desc="Processing") as pbar:
            for batch in pf.iter_batches(batch_size=BATCH_SIZE, columns=[smiles_col]):
                smiles_list = batch.column(smiles_col).to_pylist()
                # Split into chunks for parallel canonicalization
                chunk_size = max(1, len(smiles_list) // nprocs)
                chunks = [
                    smiles_list[i : i + chunk_size]
                    for i in range(0, len(smiles_list), chunk_size)
                ]
                results = pool.map(canon_batch, chunks)
                for chunk_results in results:
                    for smi in chunk_results:
                        if smi is not None:
                            bf.add(smi)
                            added += 1
                        else:
                            failed += 1
                pbar.update(len(smiles_list))

    bf.save(args.output)

    # Estimate FPR: k = (m/n) * ln(2), FPR = (1 - e^(-kn/m))^k
    m = FILTER_SIZE_BITS
    n = added
    k = (m / n) * math.log(2)
    fpr = (1 - math.exp(-k * n / m)) ** k

    print(f"\nDone! Saved to {args.output}")
    print(f"  Added:  {added:,}")
    print(f"  Failed: {failed:,}")
    print(f"  Estimated FPR: {fpr:.6e}")


if __name__ == "__main__":
    main()
