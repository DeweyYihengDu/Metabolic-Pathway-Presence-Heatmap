"""Tests for the low-level KEGG REST cache."""
from mpph.kegg import _cache_path, kegg_get


def test_cache_path_is_deterministic(tmp_path):
    assert _cache_path(tmp_path, "list/pathway") == _cache_path(tmp_path, "list/pathway")


def test_cache_path_does_not_collide_across_different_endpoints(tmp_path):
    # The readable slug alone collapses "a/b", "a//b" and "a_b" to the same
    # "a_b" (every run of non-alphanumeric characters, including a literal
    # underscore, becomes a single "_") -- these must not share a cache file.
    endpoints = ["a/b", "a//b", "a_b", "a:b", "a b"]
    paths = {_cache_path(tmp_path, e) for e in endpoints}
    assert len(paths) == len(endpoints)


def test_cache_path_readable_prefix_preserved(tmp_path):
    # The slug is still there for a human skimming the cache directory.
    assert _cache_path(tmp_path, "list/pathway").name.startswith("list_pathway_")


class _CountingSession:
    headers: dict = {}

    def __init__(self, text="response body"):
        self.calls = 0
        self.text = text

    def get(self, url, timeout=None):
        self.calls += 1
        class _Resp:
            def __init__(self, text):
                self.text = text
            def raise_for_status(self):
                pass
        return _Resp(self.text)


def test_kegg_get_uses_cache_on_second_call(tmp_path):
    session = _CountingSession()
    first = kegg_get(session, "list/pathway", tmp_path)
    second = kegg_get(session, "list/pathway", tmp_path)
    assert first == second == "response body"
    assert session.calls == 1  # second call served from disk, not re-fetched


def test_kegg_get_refresh_forces_refetch(tmp_path):
    session = _CountingSession()
    kegg_get(session, "list/pathway", tmp_path)
    kegg_get(session, "list/pathway", tmp_path, refresh=True)
    assert session.calls == 2


def test_kegg_get_write_leaves_no_temp_file_behind(tmp_path):
    session = _CountingSession()
    kegg_get(session, "list/pathway", tmp_path)
    leftovers = list(tmp_path.glob("*.tmp*"))
    assert leftovers == []
    written = list(tmp_path.glob("*.tsv"))
    assert len(written) == 1
    assert written[0].read_text(encoding="utf-8") == "response body"
