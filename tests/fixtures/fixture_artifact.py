import pytest
import main
import os

@pytest.fixture(scope="module")
def test1():
    dirpath = ("\\").join(os.path.dirname(__file__).split("\\")[:-1])
    filepath = os.path.join(dirpath, 'files\\test1.tsv')
    file = open(filepath, "r", encoding="utf8").read()
    file_lines = (file.split("\n"))
    header_line = file_lines[0]
    lines = file_lines[1:]
    yield {"file": {"header_line": header_line, "lines": lines}, "filepath": filepath}

@pytest.fixture(scope='module')
def test_client():
    flask_app = main.app
    test_client = flask_app.test_client()
    ctx = flask_app.app_context()
    ctx.push()
    yield test_client
    ctx.pop()