from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
PYTHON_HEADERS = REPOSITORY_ROOT / "include" / "cgmres" / "python"


def test_macro_continuation_backslashes_have_no_trailing_whitespace():
    offenders = []
    for header in sorted(PYTHON_HEADERS.glob("*.hpp")):
        for line_number, line in enumerate(header.read_text().splitlines(), start=1):
            stripped = line.rstrip()
            if stripped.endswith("\\") and line != stripped:
                offenders.append(f"{header.relative_to(REPOSITORY_ROOT)}:{line_number}")

    assert not offenders, (
        "Whitespace after a macro-continuation backslash: " + ", ".join(offenders)
    )
