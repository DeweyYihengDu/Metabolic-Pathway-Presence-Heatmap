"""Local, offline KO (KEGG Orthology) annotation of a protein FASTA via KOfam
HMM profiles + pyhmmer -- reimplements KofamScan's own scoring/assignment
logic (Aramaki et al. 2020, Bioinformatics), executed through pyhmmer
(Larralde et al. 2023) instead of shelling out to compiled HMMER binaries +
kofam_scan's own Ruby wrapper.

Scope: KO assignment only (no GO -- eggNOG's differently-structured
per-clade databases would be a separate, larger addition). Input must
already be gene-called/translated protein sequences (a ``.faa`` file); this
module does no ORF prediction of its own.

KOfam is real KEGG data, downloaded once (~1.5 GB compressed) via
:func:`download_kofam_db`. It is not bundled with mpph and not covered by
its MIT license -- see ``DATA_SOURCES.md``.
"""
from __future__ import annotations

import gzip
import os
import re
import shutil
import tarfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import requests

DEFAULT_KOFAM_HTTPS_BASE = "https://www.genome.jp/ftp/db/kofam/"
_KO = re.compile(r"K\d{5}")  # same pattern as mpph.userdata / mpph.kgml


@dataclass(frozen=True)
class KOListEntry:
    """One row of KOfam's ``ko_list``. A literal ``-`` in ``threshold`` or
    ``score_type`` is KEGG's own missing-value marker, parsed to ``None``."""
    ko_id: str
    threshold: float | None
    score_type: str | None
    definition: str

    @property
    def assignable(self) -> bool:
        """Some KOs have no predefined threshold (too few reference
        sequences in KEGG GENES) -- these can never be assigned via the
        significant/mapper-style output, by KOfam's own design, not a bug."""
        return self.threshold is not None


@dataclass(frozen=True)
class Assignment:
    gene_id: str
    ko_id: str
    score: float


def parse_ko_list(text: str) -> dict[str, KOListEntry]:
    """Parse KOfam's ``ko_list`` (1 header line, 12 tab-separated columns:
    knum, threshold, score_type, profile_type, F-measure, nseq, nseq_used,
    alen, mlen, eff_nseq, re/pos, definition) into ``{KO_id: KOListEntry}``.
    """
    entries: dict[str, KOListEntry] = {}
    for line in text.splitlines()[1:]:
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        ko_id, threshold_s, score_type_s = fields[0], fields[1], fields[2]
        entries[ko_id] = KOListEntry(
            ko_id=ko_id,
            threshold=None if threshold_s == "-" else float(threshold_s),
            score_type=None if score_type_s == "-" else score_type_s,
            definition=fields[11],
        )
    return entries


def load_ko_thresholds(kofam_db_dir: str | Path) -> dict[str, KOListEntry]:
    """Locate and parse ``ko_list``/``ko_list.gz`` directly under
    ``kofam_db_dir`` (either a :func:`download_kofam_db` layout or a manual
    ``gunzip ko_list.gz`` one). Raises ``FileNotFoundError`` if neither
    exists."""
    d = Path(kofam_db_dir)
    plain, gz = d / "ko_list", d / "ko_list.gz"
    if plain.exists():
        text = plain.read_text(encoding="utf-8")
    elif gz.exists():
        with gzip.open(gz, "rt", encoding="utf-8") as fh:
            text = fh.read()
    else:
        raise FileNotFoundError(
            f"ko_list not found under {d} -- run `mpph annotate --setup-db "
            f"{d}` first, or point --kofam-db at a directory containing it.")
    return parse_ko_list(text)


def discover_profiles(kofam_db_dir: str | Path) -> list[Path]:
    """Every ``*.hmm`` file under ``kofam_db_dir/profiles/`` (the layout
    both :func:`download_kofam_db` and a manual ``tar xzf profiles.tar.gz``
    produce) or directly under ``kofam_db_dir``."""
    d = Path(kofam_db_dir)
    for candidate in (d / "profiles", d):
        profiles = sorted(candidate.glob("*.hmm"))
        if profiles:
            return profiles
    raise FileNotFoundError(
        f"no *.hmm profiles found under {d} or {d / 'profiles'} -- run "
        f"`mpph annotate --setup-db {d}` first.")


def read_ko_subset(path: str | Path) -> set[str]:
    """Restrict annotation to these KO ids: one per line, ``#`` comments
    allowed. Also accepts KEGG's own ``prokaryote.hal``/``eukaryote.hal``
    (which list profile *paths* like ``profiles/K00001.hmm``, not bare ids)
    -- reusing the bare ``K#####`` scan handles both with no format-specific
    branching."""
    text = Path(path).read_text(encoding="utf-8", errors="ignore")
    kos = set(_KO.findall(text))
    if not kos:
        raise ValueError(f"no KO ids (K#####) found in {path}")
    return kos


def count_sequences(fasta_path: str | Path) -> int:
    """Number of records in a FASTA, counted without loading the sequences.

    Used to decide between pre-fetching and streaming the target database,
    so it must not itself materialise a proteome in memory.
    """
    n = 0
    with open(fasta_path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                n += 1
    return n


def is_significant_hit(entry: KOListEntry, full_sequence_score: float,
                       best_domain_score: float | None) -> bool:
    """KofamScan's own significance rule, ported from its actual source
    (``result/hit.rb``'s ``above_threshold?``) rather than its README, which
    disagrees with its own code on one point: the comparison is ``>=``, not
    "higher than". ``score_type == "domain"`` compares the best single
    domain's score; every other value (``"full"``, or undefined) compares
    the full-sequence score (``result/parser.rb``).
    """
    if entry.threshold is None:
        return False
    if entry.score_type == "domain":
        return best_domain_score is not None and best_domain_score >= entry.threshold
    return full_sequence_score >= entry.threshold


# Above this many input proteins, stream the sequence database instead of
# pre-fetching it.
#
# The number is set from a direct measurement, not from extrapolating this
# function's total memory use: a prefetched digital sequence block costs
# **~1.0 kB per protein** (measured 4.5 MB for B. subtilis' 4,237 and 45.5 MB
# for Arabidopsis' 48,265). So prefetching an entire plant proteome costs
# ~45 MB, and even the human RefSeq set (137k entries) costs ~135 MB.
#
# That is a rounding error next to this function's actual peak RSS, which is
# 2-12 GB and is dominated by the profile search itself, not by the targets
# (measured on the same input: prefetch 11.4 GB vs stream 11.6 GB -- streaming
# is, if anything, marginally worse). Streaming therefore buys nothing for any
# single organism's proteome and the threshold is placed accordingly: it bites
# only on metagenome-scale protein catalogues in the millions, where ~1 GB of
# held targets does start to matter.
PREFETCH_MAX_SEQUENCES = 1_000_000


def annotate_fasta(
    fasta_path: str | Path, kofam_db_dir: str | Path, *,
    cpus: int = 0,
    entries: dict[str, KOListEntry] | None = None,
    ko_subset: set[str] | None = None,
    prefetch: bool | None = None,
) -> list[Assignment]:
    """Search every KOfam profile in ``kofam_db_dir`` against every protein
    in ``fasta_path`` (protein FASTA, already gene-called -- no ORF
    prediction here), and return every hit for which
    :func:`is_significant_hit` is true -- i.e. this is already the final
    KofamScan-mapper-equivalent result, not raw hits.

    ``entries`` lets a caller load :func:`load_ko_thresholds` once and reuse
    it; if ``None``, this loads them itself. Needs the optional ``pyhmmer``
    dependency (``pip install mpph[annotate]``) -- imported locally so
    importing :mod:`mpph.kofam` itself never requires it.

    ``prefetch`` chooses how the target sequences are held. ``True`` loads
    them all into memory; ``False`` streams them from the file; ``None``
    (default) prefetches up to :data:`PREFETCH_MAX_SEQUENCES` proteins. Both
    paths return identical assignments -- the choice is purely about holding
    ~1 kB per protein, which no single organism's proteome makes significant
    (see :data:`PREFETCH_MAX_SEQUENCES` for the measurements). Reach for
    ``False`` on a metagenome protein catalogue, not on a genome.
    """
    import pyhmmer

    if entries is None:
        entries = load_ko_thresholds(kofam_db_dir)
    profile_paths = discover_profiles(kofam_db_dir)
    if ko_subset is not None:
        profile_paths = [p for p in profile_paths if p.stem in ko_subset]
        if not profile_paths:
            raise ValueError(
                f"--ko-subset matched none of the profiles under {kofam_db_dir}")

    alphabet = pyhmmer.easel.Alphabet.amino()
    n_sequences = count_sequences(fasta_path)
    if n_sequences == 0:
        raise ValueError(f"no protein sequences found in {fasta_path}")
    use_prefetch = (n_sequences <= PREFETCH_MAX_SEQUENCES
                    if prefetch is None else prefetch)

    def _iter_profile_hmms():
        for p in profile_paths:
            with pyhmmer.plan7.HMMFile(p) as hmm_file:
                yield from hmm_file

    assignments: list[Assignment] = []

    def _collect(hits_iter):
        # T=0: report every hit down to score 0, matching kofam_scan's own
        # `hmmsearch -T 0` invocation -- ko_list's threshold, not hmmsearch's
        # own default E-value-based gate, is meant to be the only
        # significance filter (see is_significant_hit).
        for hits in hits_iter:
            entry = entries.get(hits.query.name)
            if entry is None or not entry.assignable:
                continue
            for hit in hits:
                best_dom = hit.best_domain.score if len(hit.domains) else None
                if is_significant_hit(entry, hit.score, best_dom):
                    assignments.append(
                        Assignment(hit.name, hits.query.name, hit.score))

    if use_prefetch:
        with pyhmmer.easel.SequenceFile(fasta_path, digital=True,
                                        alphabet=alphabet) as sf:
            sequences = sf.read_block()
        _collect(pyhmmer.hmmsearch(_iter_profile_hmms(), sequences,
                                   cpus=cpus, T=0))
    else:
        with pyhmmer.easel.SequenceFile(fasta_path, digital=True,
                                        alphabet=alphabet) as sf:
            _collect(pyhmmer.hmmsearch(_iter_profile_hmms(), sf,
                                       cpus=cpus, T=0))

    assignments.sort(key=lambda a: (a.gene_id, a.ko_id))
    return assignments


def write_mapper_tsv(assignments: list[Assignment], out_path: str | Path) -> None:
    """Write ``gene<TAB>K#####`` rows, no header, sorted by (gene_id, ko_id)
    for a reproducible/diffable file -- the same shape as KofamScan
    ``--format mapper``. Already readable with zero new parsing code by
    :func:`mpph.userdata.load_user_kos`."""
    ordered = sorted(assignments, key=lambda a: (a.gene_id, a.ko_id))
    lines = [f"{a.gene_id}\t{a.ko_id}" for a in ordered]
    Path(out_path).write_text(
        "\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _extract_tar_safely(tar_path: Path, dest_dir: Path) -> None:
    with tarfile.open(tar_path) as tf:
        if hasattr(tarfile, "data_filter"):  # Python >=3.12 (or backported)
            tf.extractall(dest_dir, filter="data")
        else:  # pragma: no cover -- older Python without the safety filter
            tf.extractall(dest_dir)  # trusted KEGG source over HTTPS


def _download_one(session: requests.Session, url: str, dest_file: Path) -> None:
    """Stream ``url`` to ``dest_file`` via a temp file + atomic rename (same
    pattern as ``mpph.kegg.kegg_get``'s cache writes)."""
    tmp = dest_file.with_name(f"{dest_file.name}.tmp{os.getpid()}")
    with session.get(url, stream=True, timeout=60) as resp:
        resp.raise_for_status()
        with open(tmp, "wb") as fh:
            fh.writelines(resp.iter_content(chunk_size=1024 * 1024))
    os.replace(tmp, dest_file)


def download_kofam_db(
    dest_dir: str | Path, *,
    base_url: str = DEFAULT_KOFAM_HTTPS_BASE,
    session: requests.Session | None = None,
    overwrite: bool = False,
    progress: Callable[[str, int, int], None] | None = None,
) -> Path:
    """Download + extract ``ko_list.gz`` and ``profiles.tar.gz`` from KEGG's
    own KOfam distribution into ``dest_dir`` (created if needed). One-time
    setup, ~1.5 GB compressed. Skips re-downloading if ``dest_dir`` already
    looks populated (an existing ``ko_list``/``ko_list.gz`` and at least one
    ``profiles/*.hmm``), unless ``overwrite=True``.

    ``progress(filename, bytes_done, bytes_total)`` is called periodically
    while each file downloads, if given (``bytes_total`` is 0 if the server
    didn't send a Content-Length).
    """
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    already_populated = (
        ((dest / "ko_list").exists() or (dest / "ko_list.gz").exists())
        and any((dest / "profiles").glob("*.hmm"))
    )
    if already_populated and not overwrite:
        return dest

    sess = session or requests.Session()
    base = base_url.rstrip("/")
    for name in ("ko_list.gz", "profiles.tar.gz"):
        if progress is not None:
            progress(name, 0, 0)
        _download_one(sess, f"{base}/{name}", dest / name)

    with gzip.open(dest / "ko_list.gz", "rb") as src, \
            open(dest / "ko_list", "wb") as dst:
        shutil.copyfileobj(src, dst)
    _extract_tar_safely(dest / "profiles.tar.gz", dest)
    return dest
