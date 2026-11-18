# Contributing to chrcheck

Thank you for helping make bioinformatics reference preflight checks safer and easier to use.

## Development Setup

`chrcheck` uses only the Python 3.9+ standard library with zero runtime dependencies.

```bash
git clone https://github.com/SuzumiyaHaruhi719/chrcheck.git
cd chrcheck
python -m pip install -e .
```

## Running Tests

Run the full test suite and bytecode check locally before opening a pull request:

```bash
python -m compileall -q src tests
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Guidelines

- Keep runtime dependencies empty (standard library only) so `chrcheck` runs offline in minimal HPC and container environments.
- When adding a new reference assembly or format parser, verify contig lengths against official NCBI / UCSC / Ensembl records and add unit tests covering both passing and failing cases.
- Keep pull requests focused on a single bug fix or enhancement.
- Avoid introducing breaking changes to CLI exit codes (`0` clean/no errors, `1` compatibility error or `--fail-on warning` violation, `2` usage/parse failure) or JSON output keys.
