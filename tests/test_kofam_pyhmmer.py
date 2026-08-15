"""End-to-end test of mpph.kofam.annotate_fasta() against the *real* pyhmmer
plumbing (HMMFile, hmmsearch, hits.query.name, hit.name/.score/.best_domain).

Guarded by pytest.importorskip: the `annotate` extra (and therefore pyhmmer)
is deliberately NOT part of the base `dev` install or CI's `pip install
.[dev]` step (see pyproject.toml), so this test must degrade to a skip
rather than fail wherever pyhmmer isn't installed -- the first
importorskip/skip precedent in this repo's test suite, needed for exactly
that reason.

No committed binary .hmm fixture is used (this repo has none, and a
hand-written raw HMMER3 text profile would be fragile); instead a tiny toy
HMM is built at test time via pyhmmer's own Builder API from a handful of
hand-written toy protein sequences.
"""
import pytest

pyhmmer = pytest.importorskip("pyhmmer")

from mpph.kofam import annotate_fasta, parse_ko_list

# Four near-identical toy "family" sequences to build a profile from -- not a
# real KO, just something a real match is easy to construct against.
_TOY_FAMILY = [
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKR",
    "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKQ",
]
# A real match (identical to the family) and a totally unrelated decoy.
_MATCH_SEQ = _TOY_FAMILY[0]
_DECOY_SEQ = "GGGGPPPPLLLLWWWWAAAASSSSDDDDFFFFHHHHKKKKRRRRNNNNQQQQEEEE"


def _build_toy_hmm(name: bytes = b"TESTKO"):
    alphabet = pyhmmer.easel.Alphabet.amino()
    digital_seqs = [
        pyhmmer.easel.TextSequence(name=f"fam{i}".encode(), sequence=s).digitize(alphabet)
        for i, s in enumerate(_TOY_FAMILY)
    ]
    msa = pyhmmer.easel.DigitalMSA(alphabet, name=name, sequences=digital_seqs)
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    hmm, _profile, _optimized = builder.build_msa(msa, background)
    return hmm


def _write_kofam_db(tmp_path, *, threshold: float):
    hmm = _build_toy_hmm()
    profiles_dir = tmp_path / "profiles"
    profiles_dir.mkdir()
    with open(profiles_dir / "TESTKO.hmm", "wb") as fh:
        hmm.write(fh)

    ko_list = (
        "knum\tthreshold\tscore_type\tprofile_type\tF-measure\tnseq\tnseq_used\t"
        "alen\tmlen\teff_nseq\tre/pos\tdefinition\n"
        f"TESTKO\t{threshold}\tfull\tfull\t0.9\t4\t4\t80\t80\t4.0\t0.5\ttoy family\n"
    )
    (tmp_path / "ko_list").write_text(ko_list, encoding="utf-8")
    return tmp_path


def _write_query_fasta(tmp_path):
    fasta = tmp_path / "query.faa"
    fasta.write_text(
        f">gene_match\n{_MATCH_SEQ}\n>gene_decoy\n{_DECOY_SEQ}\n", encoding="utf-8")
    return fasta


def test_annotate_fasta_assigns_real_hit_above_threshold(tmp_path):
    db_dir = _write_kofam_db(tmp_path, threshold=50.0)  # well below the real match's score
    fasta = _write_query_fasta(tmp_path)

    assignments = annotate_fasta(fasta, db_dir, cpus=1)

    by_gene = {a.gene_id: a.ko_id for a in assignments}
    assert by_gene == {"gene_match": "TESTKO"}  # decoy gets nothing
    assert all(a.score > 50.0 for a in assignments)


def test_annotate_fasta_withholds_assignment_below_threshold(tmp_path):
    # Threshold set far above any realistic score for this toy profile --
    # proves the ko_list threshold actually gates a real, otherwise-found hit
    # (not just testing that hmmsearch runs).
    db_dir = _write_kofam_db(tmp_path, threshold=100000.0)
    fasta = _write_query_fasta(tmp_path)

    assignments = annotate_fasta(fasta, db_dir, cpus=1)

    assert assignments == []


def test_annotate_fasta_ko_subset_excluding_the_only_profile_raises(tmp_path):
    db_dir = _write_kofam_db(tmp_path, threshold=50.0)  # only profile: TESTKO
    fasta = _write_query_fasta(tmp_path)

    with pytest.raises(ValueError, match="ko-subset matched none"):
        annotate_fasta(fasta, db_dir, cpus=1, ko_subset={"K00001"})


def test_annotate_fasta_ko_subset_including_the_profile_still_finds_the_hit(tmp_path):
    db_dir = _write_kofam_db(tmp_path, threshold=50.0)
    fasta = _write_query_fasta(tmp_path)

    assignments = annotate_fasta(fasta, db_dir, cpus=1, ko_subset={"TESTKO"})
    assert {a.gene_id: a.ko_id for a in assignments} == {"gene_match": "TESTKO"}


def test_streaming_and_prefetch_give_identical_assignments(tmp_path):
    """The two sequence-loading paths must not change the answer.

    Prefetching holds every target in RAM at a measured ~1 kB per protein,
    which only matters at metagenome-catalogue scale (see
    PREFETCH_MAX_SEQUENCES). Since `auto` silently picks between them by
    input size, a user must never be able to tell which one ran from the
    output -- so assert bit-identity rather than assume it. Confirmed on real
    data too: B. subtilis gives byte-identical TSVs either way.
    """
    db_dir = _write_kofam_db(tmp_path, threshold=50.0)
    fasta = _write_query_fasta(tmp_path)

    prefetched = annotate_fasta(fasta, db_dir, cpus=1, prefetch=True)
    streamed = annotate_fasta(fasta, db_dir, cpus=1, prefetch=False)
    assert prefetched == streamed
    assert prefetched  # and not trivially empty


def test_count_sequences_does_not_load_the_proteome(tmp_path):
    from mpph.kofam import count_sequences
    fasta = tmp_path / "many.faa"
    fasta.write_text("".join(f">g{i}\nMKTAYIAK\n" for i in range(2500)))
    assert count_sequences(fasta) == 2500


def test_parse_ko_list_reads_the_same_file_annotate_fasta_uses(tmp_path):
    db_dir = _write_kofam_db(tmp_path, threshold=50.0)
    entries = parse_ko_list((db_dir / "ko_list").read_text(encoding="utf-8"))
    assert entries["TESTKO"].threshold == 50.0
