import os
import subprocess
import sys


def test():
    with open("log_test_stderr.log", "w") as log_writer:
        p = subprocess.Popen(
            [
                sys.executable,
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "libfrontend_sample.py"
                ),
                "b"
            ],
            stderr=log_writer,
            text=True,
            env={
                "LOG_LEVEL": "INFO",
                "LOG_FILE_NAME": "log_test_file.log",
                **os.environ
            }
        )
        assert p.wait() == 0
