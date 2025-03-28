from .version import __version__
from molbloom.bloom import BloomFilter, CustomFilter
import os
import molbloom.data
from importlib_resources import files
from dataclasses import dataclass
from .common import _common

_filters = {
    "zinc20": None,
    "zinc-instock": None,
    "zinc-instock-mini": None,
    "surechembl": None,
}
_descriptions = {
    "zinc20": "All ZINC20 (1,006,651,037 mols) from Oct 2021. FPR of 0.003. Requires download",
    "zinc-instock": "ZINC20 instock (9,227,726 mols). FPR of 0.0003. Requires download",
    "zinc-instock-mini": "ZINC20 instock (9,227,726 mols). FPR of 0.07. Included in package",
    "surechembl": "SureChEMBL (22,843,364 mols). FPR of 0.000025. Requires download",
}
# just put in .cache
_DEFAULT_PATH = os.path.join(os.path.expanduser("~"), ".cache", "molbloom")
_filter_urls = {
    "zinc20": "https://www.dropbox.com/s/mvn1ij9ooq5ikk9/zinc20.bloom?dl=1",
    "zinc-instock": "https://www.dropbox.com/scl/fi/vrbo5sxxr30kvg1k1m4r7/zinc-instock.bloom?rlkey=w6q8tumnkv7pqffyyq5ujzi04&st=03tcxa3i&dl=1",
    "zinc-instock-mini": None,
    "surechembl": "https://www.dropbox.com/s/f6m2wjxq42avl50/sureblcanon.bloom?dl=1",
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


def canon(smiles: str) -> str | None:
    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("To canonicalize SMILES, rdkit is required.")
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToSmiles(m, canonical=True, isomericSmiles=True, kekuleSmiles=True)


def catalogs():
    """Returns a list of available catalogs and their descriptions"""
    return _descriptions


def buy(smiles, catalog="zinc-instock", canonicalize=False, check_common=True):
    """Returns True if the SMILES is probably in the catalog, False if it is definitely not

    Parameters
    ----------
    smiles : str
        The SMILES to check
    catalog : str
        The catalog to check against. Options are 'zinc20', 'zinc-instock', 'zinc-instock-mini', 'surechembl'
    canonicalize : bool
        Whether to canonicalize the SMILES before checking
    check_common : bool
        Whether to check against a list of common reagents before checking the filter
    """
    if catalog not in _filters:
        raise ValueError(
            f"Catalog {catalog} is not available. Try one of {list(_filters.keys())}"
        )
    _load_filter(catalog)
    if canonicalize:
        smiles = canon(smiles)
        if not smiles:
            return False
    if check_common and smiles in _common:
        return True
    return smiles in _filters[catalog]
