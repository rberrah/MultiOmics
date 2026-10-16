#!/usr/bin/env python3
"""Fail-closed MOFA2 Python-provenance regression tests, no external packages."""
import copy
import hashlib
import json
import sys
from pathlib import Path

sys.path.insert(0, "scripts")
from verify_mofa_basilisk_evidence import validate

lock = json.loads(Path("multiomics-engine/environment/mofa-python-pins.json").read_text())
rows = [{"name": name, "version": version} for name, version in sorted(lock["pinnedPackages"].items())]
canonical = "".join(f"{r['name']}={r['version']}\n" for r in rows)
assert hashlib.sha256(canonical.encode()).hexdigest() == lock["fullPythonLibrarySha256"]
assert len(rows) == lock["packageCount"] == 18

baseline = {
    "schema": "multiomics.mofa-basilisk-provenance.v1",
    "source": "actual MOFA2 basiliskStart(mofa_env) Python environment, not ambient Python",
    "runMode": "verify",
    "strictLockVerified": True,
    "rPackageVersion": lock["mofaRVersion"],
    "pythonVersion": lock["pythonVersion"],
    "mofapy2Version": lock["mofapy2Target"],
    "installedPythonPackages": rows,
    "installedPackageCount": len(rows),
    "installedPythonLibrarySha256": lock["fullPythonLibrarySha256"],
    "sourceCommit": "fixture",
}
assert validate(lock, baseline, expected_commit="fixture") == lock["fullPythonLibrarySha256"]

def must_reject(label, mutate):
    x = copy.deepcopy(baseline)
    mutate(x)
    try:
        validate(lock, x, expected_commit="fixture")
    except ValueError:
        return
    raise AssertionError("Unsafe evidence accepted: " + label)

must_reject("wrong source", lambda x: x.update(source="ambient Python"))
must_reject("probe only", lambda x: x.update(runMode="probe"))
must_reject("not verified", lambda x: x.update(strictLockVerified=False))
must_reject("wrong commit", lambda x: x.update(sourceCommit="another"))
must_reject("Python version", lambda x: x.update(pythonVersion="3.12.13"))
must_reject("mofapy2 version", lambda x: x.update(mofapy2Version="0.0.0"))
must_reject("package drift", lambda x: x["installedPythonPackages"][0].update(version="0.0.0"))
must_reject("removed package", lambda x: x["installedPythonPackages"].pop())
must_reject("duplicate distribution", lambda x: x["installedPythonPackages"].append(dict(x["installedPythonPackages"][0])))
must_reject("fake checksum", lambda x: x.update(installedPythonLibrarySha256="0" * 64))
print("PASS: MOFA2 Python strict dependency lock and ten independent negative controls")
