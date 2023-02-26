from libysjs.operation import YSJSCluster


def test_run(ysjsd_test_prep):
    app, ysjsd, config = ysjsd_test_prep
    print(app.test_client().get("/ysjsd/api/v1.0/config"))
    app.test_client().post("/ysjsd/api/v1.0/stop")
