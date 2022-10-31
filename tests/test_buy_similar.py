import molbloom


def test_buy_similar():
    """Test basic build_similar interface"""
    similar_substances = molbloom.buy_similar("COC[C@H](N)[C@]1(O)CCS[C@H]1C")
    assert similar_substances.shape[0] > 1000  # last checked (11/1/2022 it is 13318)
    assert similar_substances.shape[1] == 7


def test_buy_similar_alt_databases():
    """Test basic build_similar interface"""
    similar_substances = molbloom.buy_similar(
        smiles="COC[C@H](N)[C@]1(O)CCS[C@H]1C", db="WuXi-20Q4.smi.anon"
    )
    assert similar_substances.shape[0] > 1000  # last checked (11/1/2022 it is 3172)
    assert similar_substances.shape[1] == 7


def test_buy_similar_bad_db():
    """Test basic build_similar interface"""
    import urllib

    try:
        similar_substances = molbloom.buy_similar(
            smiles="COC[C@H](N)[C@]1(O)CCS[C@H]1C", db="NOT_VALID_DATABASE"
        )
    except urllib.error.HTTPError as e:
        return
    assert "Failed to raise HTTPError" == False


def test_buy_similar_no_hits():
    """Test that a query that returns no hits has length zero"""
    similar_substances = molbloom.buy_similar("H")
    assert len(similar_substances) == 0
