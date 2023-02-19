import threading


class YSJSJob(threading.Thread):
    def __init__(self):
        super().__init__()

    def run(self):
        ...

    def get_status(self):
        ...


class YSJSCluster:
    def submit_job(self):
        ...
