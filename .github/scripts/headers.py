#!/usr/bin/env python3
"""Check paste-in C++ templates or render their PDF without editing the sources."""
from __future__ import annotations

import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
# These templates are snippets, not standalone library headers. Preserve the
# dependency order of the original workflow (its sed insertions reversed it).
DEPENDENCIES = {
    "数论.hpp": ["矩阵.hpp", "modint.hpp"],
    "图论.hpp": ["并查集.hpp", "数论.hpp"],
    "poly.hpp": ["数论.hpp"],
    "字符串.hpp": ["RMQ.hpp", "离散化.hpp", "poly.hpp"],
    "vector.hpp": ["树状数组.hpp"],
    "树.hpp": ["RMQ.hpp", "离散化.hpp"],
    "LCT.hpp": ["modint.hpp"],
    "杂项.hpp": ["dq.hpp"],
}


def read_header(path: Path) -> str:
    data = path.read_bytes()
    try:
        text = data.decode("utf-8-sig")
    except UnicodeDecodeError:
        text = data.decode("gb18030")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def headers(root: Path) -> list[Path]:
    files = sorted(root.glob("*.hpp"))
    if not files:
        raise ValueError("No root-level .hpp files found")
    for file in files:
        if file.is_symlink() or any(c in file.name for c in '\\"\n\r`'):
            raise ValueError(f"Unsupported header path: {file.name!r}")
    return files


def prepare(root: Path, target: Path) -> list[Path]:
    files = headers(root)
    names = {file.name for file in files}
    if "Header.hpp" not in names:
        raise ValueError("Missing Header.hpp")
    for file in files:
        deps = DEPENDENCIES.get(file.name, [])
        for dep in deps:
            if dep not in names:
                raise ValueError(f"{file.name}: missing dependency {dep}")
        prefix = '#pragma once\n'
        if file.name == "Header.hpp":
            prefix += '#include <bits/stdc++.h>\nusing namespace std;\n'
        else:
            prefix += '#include "Header.hpp"\n'
            prefix += ''.join(f'#include "{dep}"\n' for dep in deps)
        (target / file.name).write_text(prefix + read_header(file), encoding="utf-8")
    return files


def compile_headers(root: Path, compiler: str, jobs: int) -> int:
    if jobs < 1:
        raise ValueError("--jobs must be positive")
    executable = shutil.which(compiler)
    if not executable:
        raise ValueError(f"Compiler not found: {compiler}")
    subprocess.run([executable, "--version"], check=True)
    with tempfile.TemporaryDirectory(prefix="xcpc-check-") as directory:
        target = Path(directory)
        files = prepare(root, target)

        def check(file: Path) -> tuple[str, int, str]:
            unit = target / (file.name + ".cpp")
            unit.write_text(f'#include "{file.name}"\n', encoding="utf-8")
            try:
                result = subprocess.run(
                    [executable, "-std=c++20", "-Wall", "-Wextra", "-fsyntax-only",
                     "-fmax-errors=5", str(unit)],
                    capture_output=True, text=True, errors="replace", timeout=90,
                )
                return file.name, result.returncode, result.stdout + result.stderr
            except subprocess.TimeoutExpired:
                return file.name, 1, "Compilation exceeded 90 seconds\n"

        with ThreadPoolExecutor(max_workers=jobs) as pool:
            results = list(pool.map(check, files))
    for name, code, output in results:
        print(f"::group::{name}: {'FAIL' if code else 'PASS'}")
        print(output, end="")
        print("::endgroup::")
    failed = sum(code != 0 for _, code, _ in results)
    summary = f"{compiler}: {len(results) - failed}/{len(results)} headers passed (syntax only)."
    print(summary)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as out:
            out.write(summary + "\n")
    return int(failed != 0)


def markdown(root: Path) -> str:
    parts = [r'''---
title: "C++ Header Files"
header-includes:
  - \usepackage{xeCJK}
  - \setCJKmonofont[Scale=MatchLowercase]{Noto Sans Mono CJK SC}
  - \usepackage{fvextra}
  - \DefineVerbatimEnvironment{verbatim}{Verbatim}{breaklines=true,breakanywhere=true,breakindent=0pt}
---
''']
    for file in headers(root):
        source = read_header(file).expandtabs(4)
        fence = "`" * max(3, 1 + max((len(x) for x in re.findall(r"`+", source)), default=0))
        parts.append(f"\n## `{file.name}`\n\n{fence}cpp\n{source.rstrip()}\n{fence}\n")
    document = "".join(parts)
    # XeTeX's default trailer ID includes temporary paths; derive it from content.
    identifier = hashlib.sha256(document.encode("utf-8")).hexdigest()[:32]
    trailer = r"  - \AtBeginDvi{\special{pdf:trailerid [<%s><%s>]}}" % (identifier, identifier)
    return document.replace("header-includes:\n", "header-includes:\n" + trailer + "\n", 1)


def build(root: Path, output: Path) -> None:
    output = output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    # Do not let wall-clock PDF metadata create daily binary-only commits.
    # This normalizes timestamps, not differences between toolchain versions.
    env.setdefault("SOURCE_DATE_EPOCH", "0")
    env.setdefault("FORCE_SOURCE_DATE", "1")
    env.setdefault("TZ", "UTC")
    with tempfile.TemporaryDirectory(prefix="xcpc-pdf-") as directory:
        source = Path(directory) / "headers.md"
        source.write_text(markdown(root), encoding="utf-8")
        subprocess.run([
            "pandoc", str(source), "--standalone", "--toc", "--toc-depth=2",
            "--no-highlight", "--pdf-engine=xelatex", "--pdf-engine-opt=-no-shell-escape",
            "-V", "mainfont=DejaVu Serif", "-V", "monofont=DejaVu Sans Mono",
            "-V", "CJKmainfont=Noto Serif CJK SC", "-V", "geometry:margin=1in",
            "-o", str(output),
        ], check=True, env=env, timeout=600)
    with output.open("rb") as pdf:
        if pdf.read(5) != b"%PDF-":
            raise ValueError("Pandoc did not produce a PDF")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=["check", "build", "markdown"])
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--jobs", type=int, default=2)
    parser.add_argument("--output", type=Path, default=Path("build/headers.pdf"))
    args = parser.parse_args()
    if args.command == "check":
        return compile_headers(args.root, args.compiler, args.jobs)
    if args.command == "markdown":
        print(markdown(args.root), end="")
    else:
        build(args.root, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
