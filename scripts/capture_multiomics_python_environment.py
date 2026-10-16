#!/usr/bin/env python3
"""Record real Python runtime/dependencies; inventory is not a lockfile."""
import argparse
import hashlib
import importlib.metadata as metadata
import json
from pathlib import Path
import platform
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="multiomics-python-environment.json")
    parser.add_argument("--datasets", nargs="*", default=[],
                        help="Downloaded validation datasets to fingerprint locally")
    args = parser.parse_args()
    packages = [{"name": d.metadata.get("Name", ""), "version": d.version}
                for d in metadata.distributions()]
    packages.sort(key=lambda d: (d["name"].casefold(), d["version"]))
    files = []
    for name in args.datasets:
        p = Path(name).resolve(strict=True)
        if not p.is_file():
            parser.error(f"Not a regular dataset file: {p}")
        sha = hashlib.sha256()
        with p.open("rb") as f:
            for chunk in iter(lambda: f.read(1024 * 1024), b""):
                sha.update(chunk)
        files.append({"filename": p.name, "sha256": sha.hexdigest(),
                      "size_bytes": p.stat().st_size})
    data = {"schema_version": "1.0.0",
            "python_version": sys.version, "platform": platform.platform(),
            "packages": packages, "dataset_files": files,
            "status": "inventory_only_not_a_dependency_lock"}
    Path(args.output).write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                                 encoding="utf-8")
    print(f"Wrote actual Python inventory: {args.output}")

if __name__ == "__main__":
    main()
