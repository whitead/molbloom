from .version import __version__
from molbloom.bloom import BloomFilter, CustomFilter
import os
import molbloom.data
from importlib_resources import files
from dataclasses import dataclass

_filters = {"zinc20": None, "zinc-instock": None, "zinc-instock-mini": None}
_descriptions = {
    "zinc20": "All ZINC20 (1,006,651,037 mols) from Oct 2021. FPR of 0.003. Requires download",
    "zinc-instock": "ZINC20 instock (9,227,726 mols). FPR of 0.0003. Requires download",
    "zinc-instock-mini": "ZINC20 instock (9,227,726 mols). FPR of 0.07. Included in package",
}
# just put in .cache
_DEFAULT_PATH = os.path.join(os.path.expanduser("~"), ".cache", "molbloom")
_filter_urls = {
    "zinc20": "https://www.dropbox.com/s/mvn1ij9ooq5ikk9/zinc20.bloom?dl=1",
    "zinc-instock": "https://www.dropbox.com/s/9g5ywc2n4tzev1m/instock.bloom?dl=1",
    "zinc-instock-mini": None,
}


def _download_progress(count, block_size, total_size):
    percent = int(count * block_size * 100 / total_size)
    print("\r", end="")
    print("Downloading filter... {:d}%".format(percent), end="")
    if percent >= 100:
        print("")


def _load_big_filter(name):
    # check if it's present
    filter_path = os.path.join(_DEFAULT_PATH, f"{name}.bloom")
    if not os.path.exists(_DEFAULT_PATH):
        os.makedirs(_DEFAULT_PATH)
    if not os.path.exists(filter_path):
        print(f"Starting {name} download to cache directory {_DEFAULT_PATH}")
        import urllib.request

        urllib.request.urlretrieve(
            _filter_urls[name], filter_path, reporthook=_download_progress
        )
    if not os.path.exists(filter_path):
        raise ValueError(
            f"Filter was not able to be downloaded. Try removing the cache in {_DEFAULT_PATH}"
        )
    _filters[name] = BloomFilter(str(filter_path))


def _load_filter(name):
    global _filters
    if _filters[name] is None:
        if _filter_urls[name] is None:
            # should be in package
            ap = files(molbloom.data).joinpath(f"{name}.bloom")
            _filters[name] = BloomFilter(str(ap))
        else:
            _load_big_filter(name)


def canon(smiles):
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("To canonicalize SMILES, rdkit is required.")
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)


def catalogs():
    """Returns a list of available catalogs and their descriptions"""
    return _descriptions


def buy(smiles, catalog="zinc-instock", canonicalize=False):
    """Returns True if the SMILES is probably in the catalog, False if it is definitely not"""
    if catalog not in _filters:
        raise ValueError(
            f"Catalog {catalog} is not available. Try one of {list(_filters.keys())}"
        )
    _load_filter(catalog)
    if canonicalize:
        smiles = canon(smiles)
    return smiles in _filters[catalog]
