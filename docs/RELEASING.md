# Releasing

This repo uses `pyproject.toml` versioning.

## Checklist

1. Ensure tests + lint pass:

```bash
make test
make lint
```

2. Update the version in `pyproject.toml`.

3. (Optional) Update docs/paper outputs if needed:

```bash
make all
```

4. Tag the release:

```bash
git tag -a vX.Y.Z -m "vX.Y.Z"
git push --tags
```

## Notes

- Keep releases small and reproducible.
- Prefer documenting behavior changes in the PR description and (optionally) `docs/history/`.
