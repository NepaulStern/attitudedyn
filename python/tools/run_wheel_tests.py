#!/usr/bin/env python3
"""Run the pytest suite against the wheel installed in this interpreter.

The working directory of a cibuildwheel test is outside the checkout, so the
GSS files are found through ``ATTDYN_REPOSITORY_ROOT``.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: run_wheel_tests.py REPOSITORY_ROOT", file=sys.stderr)
        return 2
    repository_root = Path(sys.argv[1]).resolve()
    tests = repository_root / "python" / "tests"
    if not tests.is_dir():
        print(f"pytest directory not found: {tests}", file=sys.stderr)
        return 2
    os.environ["ATTDYN_REPOSITORY_ROOT"] = str(repository_root)
    return subprocess.call([sys.executable, "-m", "pytest", "-q", str(tests)])


if __name__ == "__main__":
    raise SystemExit(main())
