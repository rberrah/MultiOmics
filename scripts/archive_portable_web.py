#!/usr/bin/env python3
"""Create a deterministic, distributable ZIP of the verified browser-only release."""
import hashlib
import json
import os
from pathlib import Path
from zipfile import ZipFile, ZipInfo, ZIP_DEFLATED

root = Path("release/multiomics-web")
target = Path("release/MultiOmics-portable-web.zip")
manifest = json.loads((root / "manifest.sha256.json").read_text(encoding="utf-8"))
assert manifest["schema"] == "multiomics.portable-web.v1"
tracked = [f["path"] for f in manifest["files"]] + ["manifest.sha256.json"]
assert len(tracked) == len(set(tracked))
for relative in tracked:
    assert ".." not in Path(relative).parts and not Path(relative).is_absolute()
    path = root / relative
    assert path.is_file() and not path.is_symlink()
    if relative != "manifest.sha256.json":
        data = path.read_bytes()
        item = next(row for row in manifest["files"] if row["path"] == relative)
        assert hashlib.sha256(data).hexdigest() == item["sha256"]
        assert len(data) == item["bytes"]
target.parent.mkdir(parents=True, exist_ok=True)
with ZipFile(target, "w", compression=ZIP_DEFLATED, compresslevel=9) as out:
    for relative in sorted(tracked):
        info = ZipInfo("MultiOmics-portable-web/" + relative, date_time=(2026, 1, 1, 0, 0, 0))
        info.compress_type = ZIP_DEFLATED
        info.external_attr = (0o100644 << 16)
        out.writestr(info, (root / relative).read_bytes(), compress_type=ZIP_DEFLATED, compresslevel=9)
sha = hashlib.sha256(target.read_bytes()).hexdigest()
(target.with_suffix(target.suffix + ".sha256")).write_text(sha + "  " + target.name + "\n", encoding="utf-8")
print("PASS portable standalone ZIP: {} bytes, sha256 {}".format(target.stat().st_size, sha))
print("The ZIP contains no compiled R backend. Use the separate source checkout for Bioconductor.")
