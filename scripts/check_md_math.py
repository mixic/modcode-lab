#!/usr/bin/env python3
"""Report markdown that stops GitHub from rendering mathematics.

GitHub renders ``$...$`` and ``$$...$$`` with KaTeX only when the delimiters
survive its CommonMark parse. Three patterns break the parse silently: the
formula is shown as raw TeX text, with no error message.

* Inline ``$`` glued to the preceding character, as in ``base-$b$`` or
  ``"$a$``. The opening delimiter must follow whitespace, the start of the
  line, or ``(``.
* Inline ``$`` immediately followed by a letter, digit or underscore, as in
  ``$a$bit``. The closing delimiter must be followed by whitespace,
  punctuation, or the end of the line.
* A ``$$`` block that is indented, not separated from the surrounding text by
  blank lines, or whose body contains a line that CommonMark reads as a new
  block (a list item ``- ...`` / ``+ ...`` / ``1. ...``, a heading, a quote, a
  rule). Such a line ends the paragraph and the ``$$`` never pairs.

Usage::

    python scripts/check_md_math.py docs README.md
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

BLOCK_START = re.compile(
    r"^ {0,3}(?:[-*+](?: |$)|\d{1,9}[.)](?: |$)|#{1,6} |>|(?:-{3,}|\*{3,}|_{3,})\s*$)"
)
FENCE = re.compile(r"^\s*(`{3,}|~{3,})")

# Macros that renderers may refuse outright.  \operatorname is rejected by the
# renderer used for these notes ("The following macros are not allowed:
# operatorname"); \mathrm{...} is portable and typesets the same upright shape.
NON_PORTABLE = {r"\operatorname": r"write \mathrm{...} instead"}


def check_macros(path: Path, lines: list[str]) -> list[str]:
    """Report macros that some renderers refuse to typeset."""
    problems: list[str] = []
    for number, line in enumerate(lines, start=1):
        for macro, hint in NON_PORTABLE.items():
            if macro in line:
                problems.append(f"{path}:{number}: {macro} is not portable; {hint}")
    return problems


def mask_inline_code(line: str) -> str:
    """Blank out backtick spans, keeping the length of the line."""
    chars = list(line)
    i = 0
    while i < len(chars):
        if chars[i] != "`":
            i += 1
            continue
        run = 1
        while i + run < len(chars) and chars[i + run] == "`":
            run += 1
        j = i + run
        close = -1
        while j < len(chars):
            if chars[j] != "`":
                j += 1
                continue
            width = 1
            while j + width < len(chars) and chars[j + width] == "`":
                width += 1
            if width == run:
                close = j
                break
            j += width
        if close < 0:
            i += run
            continue
        for k in range(i, close + run):
            chars[k] = " "
        i = close + run
    return "".join(chars)


def mask_code(text: str) -> list[str]:
    """Blank out fenced blocks and inline code, keeping the line structure."""
    lines = text.splitlines()
    masked: list[str] = []
    fence_char = ""
    fence_len = 0
    for line in lines:
        match = FENCE.match(line)
        if fence_len == 0 and match:
            fence_char = match.group(1)[0]
            fence_len = len(match.group(1))
            masked.append("")
            continue
        if fence_len:
            stripped = line.strip()
            if set(stripped) == {fence_char} and len(stripped) >= fence_len:
                fence_len = 0
            masked.append("")
            continue
        masked.append(mask_inline_code(line))
    return masked


def check_inline(path: Path, number: int, line: str) -> list[str]:
    """Check one line for inline-math delimiters GitHub will not pair up."""
    problems: list[str] = []
    i = 0
    while i < len(line):
        if line[i] != "$":
            i += 1
            continue
        if line[i + 1 : i + 2] == "$":
            i += 2
            continue
        close = -1
        for j in range(i + 1, len(line)):
            if line[j] == "$" and line[j - 1] != "$" and line[j + 1 : j + 2] != "$":
                close = j
                break
        if close < 0:
            problems.append(f"{path}:{number}: opening $ has no closing $ on the line")
            break
        before = line[i - 1] if i else ""
        after = line[close + 1] if close + 1 < len(line) else ""
        if before and not before.isspace() and before != "(":
            problems.append(
                f"{path}:{number}: inline $ follows {before!r};"
                " put a space or ( before it"
            )
        if after and (after.isalnum() or after == "_"):
            problems.append(
                f"{path}:{number}: inline $ precedes {after!r};"
                " put a space or punctuation after it"
            )
        i = close + 1
    return problems


def check_file(path: Path) -> list[str]:
    """Return every GitHub math-parsing problem found in one markdown file."""
    raw = path.read_text(encoding="utf-8").splitlines()
    lines = mask_code("\n".join(raw))
    problems: list[str] = []
    in_display = False
    opened_at = 0
    for number, line in enumerate(lines, start=1):
        if line.strip() == "$$":
            indent = len(line) - len(line.lstrip())
            if indent:
                problems.append(
                    f"{path}:{number}: display $$ indented by {indent};"
                    " GitHub needs column 0"
                )
            if in_display:
                if number < len(lines) and lines[number].strip():
                    problems.append(
                        f"{path}:{number}: closing $$ needs a blank line below it"
                    )
                in_display = False
            else:
                if number > 1 and lines[number - 2].strip():
                    problems.append(
                        f"{path}:{number}: opening $$ needs a blank line above it"
                    )
                in_display = True
                opened_at = number
            continue
        if in_display:
            if line.strip() and BLOCK_START.match(line):
                problems.append(
                    f"{path}:{number}: this line ends the paragraph inside the"
                    f" $$ block opened at {opened_at}, so the math never forms"
                )
            continue
        problems.extend(check_inline(path, number, line))
    if in_display:
        problems.append(f"{path}:{opened_at}: $$ block is never closed")
    problems.extend(check_macros(path, lines))
    return problems


def find_files(targets: list[str]) -> list[Path]:
    """Expand the arguments into a sorted list of markdown files."""
    files: list[Path] = []
    for target in targets:
        path = Path(target)
        if path.is_dir():
            files.extend(sorted(path.rglob("*.md")))
        elif path.suffix == ".md" and path.is_file():
            files.append(path)
    return files


def main(argv: list[str]) -> int:
    """Check every markdown file named or found, and report the result."""
    files = find_files(argv[1:] or ["."])
    problems: list[str] = []
    for path in files:
        problems.extend(check_file(path))
    for problem in problems:
        print(problem)
    print(f"checked {len(files)} markdown files, {len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
