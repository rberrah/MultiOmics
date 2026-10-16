#!/usr/bin/env python3
"""Fail-closed numerical replay comparison; not proof of biological validity."""
import argparse
import json
import math
from pathlib import Path
import sys


def compare(left, right, path, errors, rtol, atol):
    if isinstance(left, bool) or isinstance(right, bool):
        if left is not right:
            errors.append(f"{path}: status/boolean mismatch")
    elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
        if not math.isfinite(left) or not math.isfinite(right) or not math.isclose(
            left, right, rel_tol=rtol, abs_tol=atol
        ):
            errors.append(f"{path}: numerical mismatch original={left}, replay={right}")
    elif isinstance(left, dict) and isinstance(right, dict):
        if set(left) != set(right):
            errors.append(f"{path}: missing/extra result keys {sorted(set(left)^set(right))}")
        for key in sorted(set(left) & set(right)):
            compare(left[key], right[key], f"{path}.{key}", errors, rtol, atol)
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            errors.append(f"{path}: row/list count {len(left)} vs {len(right)}")
        for index, (a, b) in enumerate(zip(left, right)):
            compare(a, b, f"{path}[{index}]", errors, rtol, atol)
    elif left != right:
        errors.append(f"{path}: original/replay values differ")


def verify(expected, actual, rtol=1e-7, atol=1e-9):
    errors = []
    methods = expected.get("methods") if isinstance(expected, dict) else None
    replayed = actual.get("methods") if isinstance(actual, dict) else None
    if not isinstance(methods, dict) or not isinstance(replayed, dict):
        return ["R backend methods missing from one or both responses"]
    if not any(isinstance(method, dict) and method.get("status") == "ok"
               for method in methods.values()):
        return ["No successful native R method to verify"]
    compare(methods, replayed, "methods", errors, rtol, atol)
    # Preprocessing is also method-defining: missing filters or changed steps
    # must be detectable, not declared an equivalent analysis.
    compare(expected.get("preprocessing"), actual.get("preprocessing"),
            "preprocessing", errors, rtol, atol)
    for field in ("analysisMode", "analysedLayers"):
        compare(expected.get(field), actual.get(field), field, errors, rtol, atol)
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path,
                        help="Extracted reproduction/ directory after R replay")
    parser.add_argument("--rtol", type=float, default=1e-7)
    parser.add_argument("--atol", type=float, default=1e-9)
    args = parser.parse_args()
    if not 0 <= args.rtol < 0.01 or not 0 <= args.atol < 0.01:
        parser.error("Tolerances must be finite and between zero and 0.01")
    expected = json.loads((args.directory/"reference_expected.json").read_text("utf-8"))
    replayed = json.loads((args.directory/"reference_replayed.json").read_text("utf-8"))
    problems = verify(expected, replayed, args.rtol, args.atol)
    if problems:
        for message in problems[:25]:
            print("FAIL:", message, file=sys.stderr)
        print(f"FAIL: {len(problems)} replay discrepancies", file=sys.stderr)
        return 1
    print("PASS: R method values, structures and preprocessing match the original "
          "reference response within declared tolerances. This does NOT validate "
          "JavaScript tables, SVG figures or biological inference.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
