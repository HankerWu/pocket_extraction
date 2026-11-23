import requests
import pytest

from pocket_extraction.data_utils import fetch_pdb


@pytest.mark.parametrize("identifier", ["../outside", "/tmp/outside", "1ab/", "1ab?", "１abc", "", None])
def test_invalid_pdb_id_preserves_files_and_skips_network(tmp_path, monkeypatch, identifier):
    output = tmp_path / "downloads"
    output.mkdir()
    previous = tmp_path / "outside.pdb"
    previous.write_bytes(b"previous structure")

    def unexpected_request(*args, **kwargs):
        pytest.fail("invalid identifiers must not reach the network")

    monkeypatch.setattr(requests, "get", unexpected_request)
    with pytest.raises(ValueError, match="PDB ID"):
        fetch_pdb(identifier, str(output), quiet=True)
    assert previous.read_bytes() == b"previous structure"
    assert not list(output.iterdir())


def response(content=b"complete structure"):
    result = requests.Response()
    result.status_code = 200
    result._content = content
    result._content_consumed = True
    return result


def test_download_does_not_remove_a_staging_file_it_could_not_open(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from pocket_extraction import data_utils

    owned_by_another_writer = tmp_path / ".1abc.pdb.fixed.tmp"
    owned_by_another_writer.write_bytes(b"another writer")
    monkeypatch.setattr(data_utils, "uuid4", lambda: SimpleNamespace(hex="fixed"))
    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: response())

    with pytest.raises(FileExistsError):
        fetch_pdb("1abc", str(tmp_path), quiet=True)
    assert owned_by_another_writer.read_bytes() == b"another writer"
    assert not (tmp_path / "1abc.pdb").exists()


def test_download_retries_timeout_and_publishes_complete_file(tmp_path, monkeypatch):
    destination = tmp_path / "1abc.pdb"
    destination.write_bytes(b"previous structure")
    attempts = []

    def get(url, *, timeout, stream):
        attempts.append((url, timeout, stream))
        assert destination.read_bytes() == b"previous structure"
        if len(attempts) < 3:
            raise requests.Timeout("read stalled")
        return response()

    monkeypatch.setattr(requests, "get", get)
    assert fetch_pdb("1abc", str(tmp_path), quiet=True, timeout=0.25) == str(destination)
    assert destination.read_bytes() == b"complete structure"
    assert attempts == [("https://files.rcsb.org/download/1abc.pdb", 0.25, True)] * 3
    assert list(tmp_path.iterdir()) == [destination]


def test_download_exhaustion_preserves_previous_file(tmp_path, monkeypatch):
    destination = tmp_path / "1abc.pdb"
    destination.write_bytes(b"previous structure")
    attempts = []

    def get(url, *, timeout, stream):
        attempts.append(timeout)
        raise requests.Timeout("connection stalled")

    monkeypatch.setattr(requests, "get", get)
    with pytest.raises(requests.Timeout):
        fetch_pdb("1abc", str(tmp_path), retry=2, quiet=True)
    assert attempts == [30.0, 30.0]
    assert destination.read_bytes() == b"previous structure"
    assert list(tmp_path.iterdir()) == [destination]


def test_interrupted_body_is_discarded_before_retry(tmp_path, monkeypatch):
    destination = tmp_path / "1abc.pdb"
    destination.write_bytes(b"previous")

    def partial_body(*, chunk_size):
        yield b"incomplete"
        raise requests.ConnectionError("body interrupted")

    interrupted = response()
    interrupted.iter_content = partial_body
    replies = iter((interrupted, response(b"valid replacement")))

    def get(*args, **kwargs):
        assert destination.read_bytes() == b"previous"
        assert list(tmp_path.iterdir()) == [destination]
        return next(replies)

    monkeypatch.setattr(requests, "get", get)
    fetch_pdb("1abc", str(tmp_path), quiet=True, retry=2)
    assert destination.read_bytes() == b"valid replacement"
    assert list(tmp_path.iterdir()) == [destination]


def test_download_publish_failure_preserves_previous_file(tmp_path, monkeypatch):
    destination = tmp_path / "1abc.pdb"
    destination.write_bytes(b"previous structure")
    monkeypatch.setattr(requests, "get", lambda *args, **kwargs: response())

    def fail_replace(source, target):
        assert source.read_bytes() == b"complete structure"
        raise OSError("publication failed")

    monkeypatch.setattr(type(destination), "replace", fail_replace)
    with pytest.raises(OSError, match="publication failed"):
        fetch_pdb("1abc", str(tmp_path), quiet=True)
    assert destination.read_bytes() == b"previous structure"
    assert list(tmp_path.iterdir()) == [destination]


@pytest.mark.parametrize("options", [
    {"retry": 0}, {"retry": -1}, {"retry": True}, {"retry": 1.5},
    {"timeout": 0}, {"timeout": float("inf")}, {"timeout": float("nan")},
])
def test_invalid_download_limits_fail_before_network(tmp_path, monkeypatch, options):
    def unexpected_request(*args, **kwargs):
        pytest.fail("invalid limits must not reach the network")

    monkeypatch.setattr(requests, "get", unexpected_request)
    with pytest.raises(ValueError):
        fetch_pdb("1abc", str(tmp_path), quiet=True, **options)
    assert not list(tmp_path.iterdir())
