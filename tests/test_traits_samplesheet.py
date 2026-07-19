"""Tests for trait panels and sample-sheet / metadata import."""
import pytest

from mpph.samplesheet import (
    import_checkm2,
    import_gtdbtk,
    load_kos_from_sheet,
    read_sample_sheet,
    sheet_metadata,
)
from mpph.traits import load_trait_panel, panel_provenance, score_trait, score_traits


def test_score_trait_all_and_any():
    steps = [{"all_of": ["K1", "K2"]}, {"any_of": ["K3", "K4"]}]
    assert score_trait(steps, {"K1", "K2", "K3"}) == 1.0
    assert score_trait(steps, {"K1", "K2"}) == 0.5
    assert score_trait(steps, {"K3"}) == 0.5


def test_optional_step_not_required():
    steps = [{"all_of": ["K1"]}, {"optional": ["K9"]}]
    assert score_trait(steps, {"K1"}) == 1.0  # optional ignored


def test_builtin_panel_loads_and_scores():
    panel = load_trait_panel("biogeochemistry")
    assert "nitrogen_fixation" in panel
    org_kos = {"diazotroph": {"K02586", "K02588", "K02591"},
               "none": {"K99999"}}
    matrix, names, cats = score_traits(org_kos, panel)
    assert matrix.loc["diazotroph", "nitrogen_fixation"] == 1.0
    assert matrix.loc["none", "nitrogen_fixation"] == 0.0
    assert cats["nitrogen_fixation"] == "Nitrogen cycle"


def test_additional_builtin_panels_valid():
    for name in ("respiration", "carbon_fixation"):
        panel = load_trait_panel(name)
        assert panel
        for _tid, defn in panel.items():
            assert defn.get("steps"), f"{name} trait missing steps"


def test_builtin_panels_have_version_and_provenance():
    for name in ("biogeochemistry", "respiration", "carbon_fixation"):
        prov = panel_provenance(name)
        assert prov["panel_version"] == "1.0.0"
        assert len(prov["panel_sha256"]) == 64
        assert prov["panel_path"].endswith(f"{name}.json")
        # panel_provenance must resolve built-in names the same way
        # load_trait_panel does, without disturbing that function's contract.
        assert load_trait_panel(name)


def test_panel_provenance_detects_content_change(tmp_path):
    panel_file = tmp_path / "custom.json"
    panel_file.write_text('{"panel_version": "1.0.0", "traits": '
                          '{"t1": {"steps": [{"all_of": ["K00001"]}]}}}')
    first = panel_provenance(panel_file)
    panel_file.write_text('{"panel_version": "1.0.1", "traits": '
                          '{"t1": {"steps": [{"all_of": ["K00002"]}]}}}')
    second = panel_provenance(panel_file)
    assert first["panel_sha256"] != second["panel_sha256"]
    assert second["panel_version"] == "1.0.1"


def test_sample_sheet_roundtrip(tmp_path):
    (tmp_path / "m1.txt").write_text("K00001\nK00002\n")
    (tmp_path / "m2.txt").write_text("K00001\n")
    sheet_file = tmp_path / "samples.tsv"
    sheet_file.write_text(
        "sample_id\tannotation_file\tinput_format\tgroup\n"
        "MAG1\tm1.txt\tko-list\tsurface\n"
        "MAG2\tm2.txt\tko-list\tdeep\n"
    )
    sheet = read_sample_sheet(sheet_file)
    kos = load_kos_from_sheet(sheet, base_dir=tmp_path)
    assert kos["MAG1"] == {"K00001", "K00002"}
    assert kos["MAG2"] == {"K00001"}
    meta = sheet_metadata(sheet)
    assert meta.loc["MAG1", "group"] == "surface"
    assert "annotation_file" not in meta.columns


def test_duplicate_sample_ids_rejected(tmp_path):
    f = tmp_path / "s.tsv"
    f.write_text("sample_id\tannotation_file\nA\tx\nA\ty\n")
    with pytest.raises(ValueError):
        read_sample_sheet(f)


def test_empty_sample_id_rejected(tmp_path):
    f = tmp_path / "s.tsv"
    f.write_text("sample_id\tannotation_file\nA\tx\n\ty\n")
    with pytest.raises(ValueError, match="empty sample_id"):
        read_sample_sheet(f)


def test_empty_annotation_file_rejected(tmp_path):
    # A blank annotation_file must not silently resolve to base_dir itself
    # (Path(base) / "" == base) and absorb every file under it.
    (tmp_path / "real_sample.txt").write_text("K00001\n")
    (tmp_path / "stray_notes.txt").write_text("unrelated K09999 mention\n")
    f = tmp_path / "s.tsv"
    f.write_text("sample_id\tannotation_file\n"
                "real_sample\treal_sample.txt\nghost_row\t\n")
    sheet = read_sample_sheet(f)
    with pytest.raises(ValueError, match="empty annotation_file"):
        load_kos_from_sheet(sheet, base_dir=tmp_path)


def test_checkm2_and_gtdbtk_import(tmp_path):
    c = tmp_path / "checkm2.tsv"
    c.write_text("Name\tCompleteness\tContamination\nMAG1\t95.4\t1.2\n")
    q = import_checkm2(c)
    assert q.loc["MAG1", "completeness"] == 95.4

    g = tmp_path / "gtdbtk.tsv"
    g.write_text("user_genome\tclassification\nMAG1\td__Bacteria;p__X\n")
    t = import_gtdbtk(g)
    assert "Bacteria" in t.loc["MAG1", "taxonomy"]
