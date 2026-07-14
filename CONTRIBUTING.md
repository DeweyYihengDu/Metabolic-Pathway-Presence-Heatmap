# Contributing

Thanks for your interest in improving MPPH.

## Development setup

```bash
git clone https://github.com/DeweyYihengDu/Metabolic-Pathway-Presence-Heatmap.git
cd Metabolic-Pathway-Presence-Heatmap
pip install -e ".[dev]"
```

## Before opening a pull request

- `ruff check mpph tests` — lint must pass.
- `pytest -q` — all tests must pass. Unit tests and the offline integration
  test (`tests/test_integration.py`) run entirely from cache fixtures and must
  not require network access.
- Add or update tests for any behaviour change, especially in the module
  completeness scorer (`mpph/modules.py`) and the matrix/filtering logic.
- Keep the CLI, `--version`, and `pyproject.toml` version in sync.

## Scope

MPPH is a KEGG functional-profile visualiser. Please keep changes focused;
larger features (reference-tree comparison, bootstrap support, interactive
reports) are welcome but are best discussed in an issue first.
