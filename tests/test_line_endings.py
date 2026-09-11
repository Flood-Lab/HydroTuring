"""Everything Git treats as text is stored with LF line endings.

.gitattributes normalizes line endings on commit, but only in a checkout that
already has it, and under text=auto Git does not convert a file whose copy in
the index already contains CRLF. A branch cut before the attribute, or a tool
that ignores it, can still commit a file as CRLF, and the pull request then
shows every line of that file as changed. #32 arrived that way.

The check reads the index, not the working tree, so a Windows checkout passes
whatever its editor writes to disk as long as what it commits is LF. ht.cmd
needs no exception for the same reason: it is CRLF on disk and LF in Git.
"""

from __future__ import annotations

import shutil
import subprocess

import pytest

from hydroturing.spec import REPO_ROOT


@pytest.fixture(scope="module")
def tracked() -> list[tuple[str, str, str]]:
    """(line endings in the index, gitattributes, path) for every tracked file."""
    if not (REPO_ROOT / ".git").exists() or shutil.which("git") is None:
        pytest.skip("not a git checkout, or git is not on PATH")
    out = subprocess.run(
        ["git", "ls-files", "--eol", "-z"],
        cwd=REPO_ROOT, capture_output=True, check=True, encoding="utf-8",
    ).stdout
    entries = []
    for record in filter(None, out.split("\0")):
        info, _, path = record.partition("\t")
        index, _, *attr = info.split()
        entries.append((index.removeprefix("i/"), " ".join(attr).removeprefix("attr/"), path))
    return entries


def test_text_is_committed_with_lf(tracked):
    crlf = [
        f"{path} ({index})"
        for index, attr, path in tracked
        if index in ("crlf", "mixed") and attr != "-text"
    ]
    assert not crlf, (
        f"committed with CRLF line endings: {', '.join(crlf)}. Merge main into the branch, "
        "run `git add --renormalize .` and commit; the diff then shows only the lines that changed"
    )


def test_windows_launcher_is_still_checked_out_with_crlf(tracked):
    attrs = {path: attr for _, attr, path in tracked}
    assert attrs["ht.cmd"] == "text eol=crlf", (
        f"ht.cmd resolves to {attrs['ht.cmd']!r}: `ht.cmd text eol=crlf` is missing from "
        ".gitattributes or a later line overrides it, so Windows would check it out with LF"
    )
