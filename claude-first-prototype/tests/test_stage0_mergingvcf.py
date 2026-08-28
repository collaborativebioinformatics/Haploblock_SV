"""Offline tests for Stage 0's bcftools + truvari workflow."""

import json
import os
import shutil
import stat
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "pipeline" / "stage_0_mergingvcf.py"


def make_fake_tools(tmp_path: Path) -> tuple[Path, Path]:
    bcftools = tmp_path / "bcftools"
    bcftools.write_text(
        """#!/usr/bin/env python3
import json, os, pathlib, shutil, sys
args = sys.argv[1:]
log = pathlib.Path(os.environ["FAKE_TOOL_LOG"])
with log.open("a") as handle:
    handle.write(json.dumps({"tool": "bcftools", "args": args}) + "\\n")
if args[0] == "query" and args[1] == "-l":
    for line in pathlib.Path(args[-1]).read_text().splitlines():
        if line.startswith("#CHROM"):
            print(line.split("\\t")[-1])
            break
elif args[0] == "view" and args[1] == "-h":
    print(pathlib.Path(args[-1]).read_text(), end="")
elif args[0] in {"sort", "annotate"}:
    out = pathlib.Path(args[args.index("-o") + 1])
    shutil.copyfile(args[-1], out)
elif args[0] == "merge":
    out = pathlib.Path(args[args.index("-o") + 1])
    out.write_text("merged\\n")
elif args[0] == "index":
    if "-n" in args:
        print("0")
    else:
        pathlib.Path(args[-1] + ".tbi").write_text("index\\n")
else:
    raise SystemExit(2)
"""
    )
    truvari = tmp_path / "truvari"
    truvari.write_text(
        """#!/usr/bin/env python3
import json, os, pathlib, sys
args = sys.argv[1:]
log = pathlib.Path(os.environ["FAKE_TOOL_LOG"])
with log.open("a") as handle:
    handle.write(json.dumps({"tool": "truvari", "args": args}) + "\\n")
out = pathlib.Path(args[args.index("-o") + 1])
out.write_text("collapsed\\n")
removed = pathlib.Path(args[args.index("-c") + 1])
removed.write_text("removed\\n")
"""
    )
    for tool in (bcftools, truvari):
        tool.chmod(tool.stat().st_mode | stat.S_IXUSR)
    return bcftools, truvari


def test_stage0_merges_and_collapses_manifest(tmp_path):
    vcf_a = tmp_path / "a.vcf"
    vcf_b = tmp_path / "b.vcf"
    vcf_a.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_A\n"
    )
    vcf_b.write_text(
        "##fileformat=VCFv4.3\n"
        "##contig=<ID=1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE_B\n"
    )
    manifest = tmp_path / "vcfs.txt"
    manifest.write_text(f"# cohort\n{vcf_a}\n\n{vcf_b}\n")
    bcftools, truvari = make_fake_tools(tmp_path)
    out_dir = tmp_path / "out"
    tool_log = tmp_path / "tool-log.jsonl"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--vcf-list",
            str(manifest),
            "--out-dir",
            str(out_dir),
            "--bcftools",
            str(bcftools),
            "--truvari",
            str(truvari),
            "--truvari-arg=--refdist",
            "--truvari-arg",
            "500",
        ],
        capture_output=True,
        text=True,
        env={
            **os.environ,
            "PATH": f"{tmp_path}:{os.environ['PATH']}",
            "FAKE_TOOL_LOG": str(tool_log),
        },
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert (out_dir / "merged.exact.vcf.gz").read_text() == "merged\n"
    assert (out_dir / "merged.collapsed.vcf.gz").read_text() == "collapsed\n"
    assert (out_dir / "merged.collapsed_removed.vcf.gz").read_text() == "removed\n"
    assert (out_dir / "merged.exact.vcf.gz.tbi").exists()
    assert (out_dir / "merged.collapsed.vcf.gz.tbi").exists()
    assert (out_dir / "merged.collapsed_removed.vcf.gz.tbi").exists()
    run = json.loads((out_dir / "run.json").read_text())
    assert run["stage"] == 0
    assert run["truvari_args"] == ["--refdist", "500"]
    assert run["truvari_removed_output"].endswith("merged.collapsed_removed.vcf.gz")

    commands = [json.loads(line) for line in tool_log.read_text().splitlines()]
    merge = next(command for command in commands if command["args"][0] == "merge")
    assert merge["args"][1:3] == ["-m", "id"]
    truvari_index = next(
        index
        for index, command in enumerate(commands)
        if command["tool"] == "truvari"
    )
    output_steps = commands[truvari_index + 1 :]
    assert [command["args"][0] for command in output_steps[:4]] == [
        "sort",
        "index",
        "sort",
        "index",
    ]
    assert all(
        "--threads" not in command["args"]
        for command in output_steps
        if command["args"][0] == "sort"
    )


@pytest.mark.parametrize(
    ("manifest_text", "message"),
    [("", "contains no input VCFs"), ("{path}\n{path}\n", "duplicate input paths")],
)
def test_manifest_validation(tmp_path, manifest_text, message):
    vcf = tmp_path / "one.vcf"
    vcf.write_text("##fileformat=VCFv4.3\n")
    manifest = tmp_path / "vcfs.txt"
    manifest.write_text(manifest_text.format(path=vcf))
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--vcf-list", str(manifest)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert message in result.stderr
