"""Verify every publication file against the portable SHA256 inventory."""
import hashlib
import json
from pathlib import Path
import sys


def main(root):
    root = root.resolve()
    inventory = json.loads((root / "FILE_MANIFEST_SHA256.json").read_text(encoding="utf-8-sig"))
    expected = set()
    for record in inventory["files"]:
        path = (root / record["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Inventory path escapes the publication directory.")
        assert path.is_file(), record["path"]
        assert path.stat().st_size == record["bytes"], record["path"]
        with path.open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        assert actual == record["sha256"], record["path"]
        expected.add(record["path"])
    actual_paths = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    assert actual_paths == expected | {"FILE_MANIFEST_SHA256.json"}, "Unlisted or missing publication files."
    provenance = json.loads((root / "data_provenance.json").read_text(encoding="utf-8-sig"))
    assert not provenance["missing_raw_data"]
    for group in provenance["groups"]:
        for key in ["direct_sources", "raw_sources", "code_sources"]:
            assert set(group[key]) <= expected, (group["id"], key)
    print(f"Verified {len(expected)} file hashes, complete inventory coverage, and all provenance links.")


if __name__ == "__main__":
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent)
