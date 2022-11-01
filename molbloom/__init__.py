from .version import __version__
from molbloom.bloom import BloomFilter, CustomFilter
import os
import molbloom.data
from importlib_resources import files
from dataclasses import dataclass

_filters = {"zinc20": None, "zinc-instock": None, "zinc-instock-mini": None}
_descriptions = {
    "zinc20": "All ZINC20 (1,006,651,037 mols) from Oct 2021. FPR of 0.003. Requires download",
    "instock": "ZINC20 instock (9,227,726 mols). FPR of 0.0003. Requires download",
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


@dataclass
class SmallWorldHit:
    """Small World Similarity Search Hit Data"""

    #: Hit smiles
    smiles: str
    #: Hit compound id
    compound_id: str
    #: Hit graph edit distance to query
    dist: int
    #: Hit extended connectivity fingerprint (radius = 2) to query
    ecfp4: float
    #: Hit daylight fingerprint distance
    daylight: float
    #: Hit Maximum Common Edge Subgraph
    mces: int


def buy_similar(
    smiles,
    db="REAL-Database-22Q1.smi.anon",
    small_world_args={
        "dist": 4,
        "sdist": 12,
        "tdn": 6,
        "tup": 6,
        "rdn": 6,
        "rup": 2,
        "ldn": 2,
        "lup": 2,
        "maj": 6,
        "min": 6,
        "sub": 6,
        "scores": "Atom%20Alignment,ECFP4,Daylight",
    },
    n_retries=5,
    verbose=False,
):

    import urllib
    import urllib.parse
    import urllib.request
    import time

    try:
        sw_server_path = "https://sw.docking.org/search/"
        args = (("smi", smiles), ("db", db)) + tuple(small_world_args.items())
        query_url = f"{sw_server_path}submit?{urllib.parse.urlencode(args)}"
    except:
        raise Exception(
            f"Failed to construct sw.docking.org query url for smiles '{smiles}'"
        )

    if verbose:
        print(f"Querying ZINC Small World with url: {query_url}")

    lines = None
    http_status = None
    for attempt_i in range(n_retries):
        try:
            with urllib.request.urlopen(query_url) as response:
                http_status = response.status
                if http_status == 200:
                    lines = [line.decode("utf-8")[:-1] for line in response.readlines()]
                    break
        except urllib.error.HTTPError as e:
            print(
                f"ERROR: Failed to query https://sw.docking.org with smiles '{smiles}'"
            )
            print(f"ERROR: Query URL: {query_url}")
            print(f"ERROR: {e}")
            raise e

        if verbose:
            print(
                f"HTTP request attempt {attempt_i} failed with status: {http_status} ..."
            )
        time.sleep(1)

    if lines is None:
        raise Exception(f"Failed to query sw.docking.org, HTTPS status: {http_status}")

    sw_status = None
    hlid = None
    for line in lines:
        if line == "":
            continue
        line = line.replace("data:{", "").replace("}\n", "")
        line = line.split(",")

        for key_value in line:
            if '"status":' in key_value:
                sw_status = key_value.replace('"status":', "").replace('"', "")

            if "hlid" in key_value:
                hlid = key_value.replace('"hlid":', "")

    if sw_status is None:
        response_str = "\n".join(lines)
        raise Exception(f"Unexpected result from SmallWorld:\n{response_str}")
    elif sw_status == "MISS":
        if verbose:
            print(f"No hits found for smiles {smiles}")
        return []
    elif sw_status != "END":
        raise Exception(f"Unexpected status from SmallWorld '{sw_status}'")

    try:
        hlid = int(hlid)
    except:
        raise Exception(
            f"Expected small world query id to be an integer, instead it was {hlid}"
        )

    results_args = (
        "&".join(
            [
                f"hlid={hlid}",
                "order[0][column]=0",
                "columns[0][name]=alignment",
                "order[0][dir]=asc",
                "columns[1][name]=dist",
                "columns[1][search][value]=0-12",
                "columns[2][name]=ecfp4",
                "columns[3][name]=daylight",
                "columns[5][name]=mces",
            ]
        )
        .replace("[", "%5B")
        .replace("]", "%5D")
    )
    results_url = f"{sw_server_path}export?{results_args}"

    if verbose:
        print(f"Getting results from ZINC Small World with url: {results_url}")

    http_status = None
    hits = []
    for attempt_i in range(n_retries):
        try:
            with urllib.request.urlopen(results_url) as response:
                if response is None:
                    time.sleep(1)
                    continue
                http_status = response.status
                if http_status == 200:
                    next(response)
                    for line in response.readlines():
                        line = line.decode("utf-8")[:-1].split("\t")
                        smiles, compound_id = line[0].split(" ")
                        hits.append(
                            SmallWorldHit(
                                smiles=smiles,
                                compound_id=compound_id,
                                dist=int(line[1]),
                                ecfp4=float(line[2]),
                                daylight=float(line[3]),
                                #
                                mces=int(line[5]),
                            )
                        )
                    break
        except urllib.error.HTTPError as e:
            print(
                f"ERROR: Failed to retrieve results from https://sw.docking.org with smiles '{smiles}'"
            )
            print(f"ERROR: Results url {results_url}")
            print(f"ERROR: {e}")
            raise e

        if verbose:
            print(
                f"HTTP request attempt {attempt_i} failed with status: {http_status} ..."
            )
        time.sleep(1)

    if http_status != 200:
        raise Exception(f"Failed to get results with HTTP status: {http_status}")

    return hits
