# Tool Scripts

Scripts for building bloom filters from compound databases.

## build_surechembl.py

Builds the SureChEMBL bloom filter from a parquet export of the SureChEMBL database.

**Requirements:** `pyarrow`, `tqdm`, `rdkit`, `molbloom` (installed in editable mode)

```bash
# Install deps
uv pip install pyarrow tqdm

# Download compounds.parquet from SureChEMBL (~3.9GB)
# Then build the filter:
python tool/build_surechembl.py compounds.parquet --output sure-chembl.bloom
```

Options:
- `--output` — output bloom file path (default: `sure-chembl.bloom`)
- `--processes` — number of worker processes (default: cpu_count - 2)

The script auto-detects the SMILES column, canonicalizes all SMILES in parallel, and writes a 200MB bloom filter. Invalid SMILES are skipped and counted.

After building, upload the `.bloom` file to Dropbox and update the URL in `molbloom/__init__.py`.

## canonicalize.py

Canonicalizes `.smi` files using RDKit via multiprocessing.

```bash
python tool/canonicalize.py *.smi
```

## build_zinc.sh / build_instock.sh

Shell scripts for building the ZINC20 and ZINC instock bloom filters using the C tool (`main.c`).
