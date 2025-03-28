import molbloom
import os


def _randomize_smiles(smi):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smi)
    canonical = Chem.MolToSmiles(mol, canonical=True)
    cur = canonical
    for i in range(10):
        cur = Chem.MolToSmiles(mol, canonical=False, doRandom=True)
        if cur != canonical:
            return cur
    return None


def test_version():
    assert molbloom.__version__


def test_example():
    assert molbloom.buy(
        "Nc1ccn([C@H]2OC(CO)=C[C@@H]2O)c(=O)n1",
        catalog="zinc-instock-mini",
        canonicalize=True,
    )
    assert not molbloom.buy(
        "c1c(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)", catalog="zinc20", canonicalize=True
    )
    assert not molbloom.buy("ZZZ", catalog="zinc-instock-mini")


def test_common():
    assert molbloom.buy("Br")
    assert not molbloom.buy("Br", catalog="zinc20", check_common=False)


def test_alot():
    import timeit

    r = timeit.timeit(
        lambda: molbloom.buy("Nc1ccn([C@H]2OC(CO)=C[C@@H]2O)c(=O)n1"), number=100000
    )
    print("Timing per call: {:.0f}ns".format(r / 100000 * 1e9))


def test_fpr():
    """See how many molecules are false positives by isomerizing the SMILES"""
    count = 0
    fp = 0
    with open(os.path.join(os.path.dirname(__file__), "some_mols_instock.txt")) as f:
        for line in f:
            s = line.split()[0]
            assert molbloom.buy(s, catalog="zinc-instock-mini", canonicalize=True)
            for i in range(10):
                rs = _randomize_smiles(s)
                if rs is not None:
                    if molbloom.buy(rs, catalog="zinc-instock-mini", canonicalize=True):
                        fp += 1
                    count += 1
    assert count > 1000
    print("False positive rate for instock mini {:f} (N={})".format(fp / count, count))

    count = 0
    fp = 0
    with open(os.path.join(os.path.dirname(__file__), "some_mols_zinc.txt")) as f:
        for line in f:
            s = line.split()[0]
            assert molbloom.buy(s, catalog="zinc20")
            for i in range(10):
                rs = _randomize_smiles(s)
                if rs is not None:
                    if molbloom.buy(rs, catalog="zinc20"):
                        fp += 1
                    count += 1
    assert count > 1000
    print("False positive rate for zinc all is {:f} (N={})".format(fp / count, count))


def test_build_custom():
    """Build a custom filter"""
    bf = molbloom.CustomFilter(100, 1000, "test")
    bf.add("CCCO")
    # canonicalize
    s = molbloom.canon("CCCOC")
    bf.add(s)
    assert "CCCO" in bf
    assert "CCCOO" not in bf
    bf.save("test.bloom")

    bf = molbloom.BloomFilter("test.bloom")
    assert "CCCO" in bf
    assert "CCCOO" not in bf


def test_catalogs():
    d = molbloom.catalogs()
    assert "zinc20" in d


def test_all_catalogs():
    for catalog in molbloom.catalogs():
        molbloom.buy("CCCOC", catalog=catalog)
