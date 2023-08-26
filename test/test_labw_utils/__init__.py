import os

if os.name == "nt":
    NULL_PATH = "NUL"
else:
    NULL_PATH = "/dev/null"
