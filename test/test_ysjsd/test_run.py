def test_run(ysjsd_test_prep):
    _, _, client = ysjsd_test_prep
    resp = client.get("/")
    assert resp.status_code == 200
