import os
import subprocess
import sys


def test():
    with open("log_test_stderr.log", "w") as log_writer:
        cmd = [
            sys.executable,
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                "libfrontend_sample.py"
            ),
            "b"
        ]
        log_writer.write(f"EXEC {' '.join(cmd)}\n")
        log_writer.flush()
        p = subprocess.Popen(
            cmd,
            stderr=log_writer,
            text=True,
            env={
                "LOG_LEVEL": "INFO",
                "LOG_FILE_NAME": "log_test_file.log",
                **os.environ
            }
        )
        assert p.wait() == 0
