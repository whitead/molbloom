from .version import __version__
from molbloom.bloom import BloomFilter
import os

filter = None
# just put in .cache
DEFAULT_PATH = os.path.join(os.path.expanduser("~"), ".cache", "molbloom")
FILTER_URL = "https://www.dropbox.com/s/mvn1ij9ooq5ikk9/zinc20.bloom?dl=1"


def _download_progress(count, block_size, total_size):
    percent = int(count * block_size * 100 / total_size)
    print("\r", end="")
    print("Downloading filter... {:d}%".format(percent), end="")
    if percent >= 100:
        print("")


def _load_filter():
    global filter
    if filter is None:
        # check if it's present
        filter_path = os.path.join(DEFAULT_PATH, "zinc20.bloom")
        if not os.path.exists(DEFAULT_PATH):
            os.makedirs(DEFAULT_PATH)
            if not os.path.exists(filter_path):
                print("Downloading filter...")
                import urllib.request

                urllib.request.urlretrieve(
                    FILTER_URL, filter_path, reporthook=_download_progress
                )

        filter = BloomFilter(str(filter_path))


def buy(smiles, canonicalize=False):
    _load_filter()
    if canonicalize:
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("To canonicalize SMILES, rdkit is required.")
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonicalize=True)
    return smiles in filter
