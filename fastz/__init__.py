from .version import __version__
from importlib_resources import files
import fastz.filters
from fastz.bloom import BloomFilter

filter = None


def _load_filter():
    global filter
    if filter is None:
        ap = files(fastz.filters).joinpath("zinc20.bloom")
        filter = BloomFilter(str(ap))


def buy(smiles, canonicalize=False):
    _load_filter()
    if canonicalize:
        try:
            from rdkit import Chem
        except ImportError:
            raise ImportError("To canonicalize SMILES, rdkit is required.")
        smiles = Chem.MolToSmiles(
            Chem.MolFromSmiles(smiles), canonicalize=True)
    return smiles in filter
