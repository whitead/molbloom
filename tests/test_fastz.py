import fastz


def test_version():
    assert fastz.__version__


def test_example():
    assert fastz.buy("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)N")
    assert not fastz.buy("c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)He")
    assert not fastz.buy("ZZZ")


def test_alot():
    import timeit
    r = timeit.timeit(lambda: fastz.buy(
        "c1ccc(c(c1)C(=O)OCC[C@@H]2CCCC[NH2+]2)N"), number=100000)
    print('Timing per call: {:.0f}ns'.format(r / 100000 * 1e9))
