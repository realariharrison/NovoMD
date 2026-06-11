"""Tests for the `novomd` CLI presentation (no RDKit needed)."""

import pytest

from novomd.__about__ import __version__
from novomd.cli import main


def test_bare_invocation_shows_panel(capsys):
    """A bare `novomd` greets with the panel and exits 0 (not an argparse error)."""
    rc = main([])
    out = capsys.readouterr().out

    assert rc == 0
    assert "novomd" in out
    assert "MOLECULAR PROPERTY CALCULATOR" in out
    assert "props" in out
    assert "explain" in out
    assert "report" in out
    assert "batch" in out
    assert "novomcp.com" in out
    assert __version__ in out


def test_no_ansi_when_not_a_tty(capsys):
    """Captured (non-TTY) output carries no ANSI escapes, so pipes stay clean."""
    main([])
    out = capsys.readouterr().out
    assert "\x1b[" not in out


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        main(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_unknown_command_still_errors(capsys):
    """An invalid subcommand remains a hard argparse error (exit 2)."""
    with pytest.raises(SystemExit) as exc:
        main(["nonsense"])
    assert exc.value.code == 2
