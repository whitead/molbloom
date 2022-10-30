from .version import __version__
from molbloom.bloom import BloomFilter
import os
import molbloom.data
from importlib_resources import files

filters = {"zinc20": None, "instock": None}
# just put in .cache
DEFAULT_PATH = os.path.join(os.path.expanduser("~"), ".cache", "molbloom")
filter_urls = {
    "zinc20": "https://www.dropbox.com/s/mvn1ij9ooq5ikk9/zinc20.bloom?dl=1",
    "instock": None,
}


def _download_progress(count, block_size, total_size):
    percent = int(count * block_size * 100 / total_size)
    print("\r", end="")
    print("Downloading filter... {:d}%".format(percent), end="")
    if percent >= 100:
        print("")


def _load_big_filter(name):
    # check if it's present
    filter_path = os.path.join(DEFAULT_PATH, f"{name}.bloom")
    if not os.path.exists(DEFAULT_PATH):
        os.makedirs(DEFAULT_PATH)
    if not os.path.exists(filter_path):
        print(f"Starting filter download to cache directory {DEFAULT_PATH}")
        import urllib.request

        urllib.request.urlretrieve(
            filter_urls[name], filter_path, reporthook=_download_progress
        )
    if not os.path.exists(filter_path):
        raise ValueError(
            f"Filter was not able to be downloaded. Try removing the cache in {DEFAULT_PATH}"
        )
    filters[name] = BloomFilter(str(filter_path))


def _load_filter(name):
    global filters
    if filters[name] is None:
        if filter_urls[name] is None:
            # should be in package
            ap = files(molbloom.data).joinpath(f"{name}.bloom")
            filters[name] = BloomFilter(str(ap))
        else:
            _load_big_filter(name)


def buy(smiles, instock=True, canonicalize=False):
    if instock:
        name = "instock"
    else:
        name = "zinc20"
    _load_filter(name)
    if canonicalize:
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("To canonicalize SMILES, rdkit is required.")
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=True)
    return smiles in filters[name]
