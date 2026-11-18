"""Replay the independent exact and numerical checks and compare the results."""
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys


CHECKS = Path(__file__).resolve().parent


def run_checker(filename):
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)
    result = subprocess.run(
        [sys.executable, "-B", filename],
        cwd=CHECKS,
        env=env,
        check=True,
        stdout=subprocess.PIPE,
        text=True,
    )
    return result.stdout.splitlines()


def main():
    reference = json.loads((CHECKS / "reference-results.json").read_text())
    formal = run_checker("formal_checks.py")
    if formal != reference["formal"]:
        raise SystemExit("The exact cyclotomic results differ from the reference.")
    print("All eight formal determinants and the exact norm constants pass.", flush=True)

    signed = run_checker("signed_checks.py")
    if len(signed) != len(reference["signed"]):
        raise SystemExit("The numerical checker omitted a reference sequence.")

    # Roundoff can vary with the multiple-precision arithmetic backend.
    separator = ", maximum relative discrepancy "
    for line, expected in zip(signed, reference["signed"]):
        case, error = line.split(separator)
        expected_case, _ = expected.split(separator)
        if case != expected_case or not Decimal(0) <= Decimal(error) < Decimal("1e-88"):
            raise SystemExit(f"The signed matrix comparison failed: {line}")
        print(line)
    print("The two embedding evaluations pass for all twelve numerical cases.")


if __name__ == "__main__":
    main()
