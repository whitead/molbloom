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
    assert molbloom.buy("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)N")
    assert not molbloom.buy("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)He")
    assert not molbloom.buy("ZZZ")


def test_alot():
    import timeit

    r = timeit.timeit(
        lambda: molbloom.buy("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)N"), number=100000
    )
    print("Timing per call: {:.0f}ns".format(r / 100000 * 1e9))


def test_fpr():
    """See how many molecules are false positives by isomerizing the SMILES"""
    count = 0
    fp = 0
    with open(os.path.join(os.path.dirname(__file__), "some_mols.txt")) as f:
        for line in f:
            s = line.split()[0]
            assert molbloom.buy(s)
            for i in range(10):
                rs = _randomize_smiles(s)
                if rs is not None:
                    if molbloom.buy(rs):
                        fp += 1
                    count += 1
    assert count > 1000
    print("False positive rate is {:f} (N={})".format(fp / count, count))
