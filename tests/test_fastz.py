import fastz


def test_version():
    assert fastz.__version__


def test_example():
    assert fastz.purchasable("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)N")
    assert not fastz.purchasable("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)He")
    assert not fastz.purchasable("ZZZ")
