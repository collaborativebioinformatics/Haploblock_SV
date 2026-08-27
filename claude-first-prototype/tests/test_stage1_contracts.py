"""Contract tests for Stage 1's downstream tables."""

import csv
import gzip
import sys
from pathlib import Path

import yaml

PIPELINE_DIR = Path(__file__).resolve().parents[1] / "pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

import stage1_cluster_aware
import match_svs_to_clusters


def write_test_vcf(path: Path) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGM00001\tSAMPLE2",
        "1\t101\tsv1\tN\t<DEL>\t.\tPASS\tEND=200;SVTYPE=DEL;SVLEN=-100\tGT\t0/1\t1/1",
        "1\t301\tsv2\tN\t<INS>\t.\tPASS\tEND=302;SVTYPE=INS;SVLEN=25\tGT\t0/0\t./.",
    ]
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage1_publishes_downstream_contract(tmp_path: Path) -> None:
    vcf_path = tmp_path / "cohort.vcf.gz"
    write_test_vcf(vcf_path)

    metadata_path = tmp_path / "sample_metadata.tsv"
    metadata_path.write_text(
        "sample_id\tpopulation\tsuperpopulation\n"
        "EXTRA\tOTHER\tOTHER\n"
        "SAMPLE2\tPOP2\tSUPER2\n"
        "GM00001\tPOP1\tSUPER1\n"
    )

    cluster_dir = tmp_path / "clusters" / "chr1"
    cluster_dir.mkdir(parents=True)
    (cluster_dir / "chr1_0-1000_cluster.tsv").write_text(
        "cluster_a\tGM00001_chr1_region_0-1000_hap0\n"
        "cluster_b\tGM00001_chr1_region_0-1000_hap1\n"
        "cluster_a\tSAMPLE2_chr1_region_0-1000_hap0\n"
        "cluster_b\tSAMPLE2_chr1_region_0-1000_hap1\n"
    )

    out_dir = tmp_path / "stage1_output"
    stage1_cluster_aware.main(
        [
            "--vcf", str(vcf_path),
            "--sample-metadata", str(metadata_path),
            "--cluster-root", str(tmp_path / "clusters"),
            "--chroms", "chr1",
            "--chrom-workers", "1",
            "--download-workers", "1",
            "--out-dir", str(out_dir),
        ]
    )

    config = yaml.safe_load((out_dir / "config.yaml").read_text())
    assert set(config["paths"]) == {
        "vcf",
        "samples",
        "sample_metadata",
        "sv_genotypes",
        "sv_to_clusters",
        "haploblocks",
        "cluster_memberships",
        "sv_block_summary",
        "debug_and_qc",
    }
    for key in ("sv_genotypes", "sv_to_clusters", "haploblocks", "cluster_memberships", "sv_block_summary"):
        assert Path(config["paths"][key]["chr1"]).exists()

    samples = read_tsv(out_dir / "samples.tsv")
    assert samples == [
        {"sample_id": "NA00001", "original_sample_id": "GM00001"},
        {"sample_id": "SAMPLE2", "original_sample_id": "SAMPLE2"},
    ]

    normalized_metadata = read_tsv(out_dir / "sample_metadata.tsv")
    assert normalized_metadata == [
        {"sample_id": "NA00001", "population": "POP1", "superpopulation": "SUPER1"},
        {"sample_id": "SAMPLE2", "population": "POP2", "superpopulation": "SUPER2"},
    ]

    genotypes = read_tsv(out_dir / "sv_genotypes.chr1.tsv")
    assert [row["sv_id"] for row in genotypes] == ["sv1", "sv2"]
    assert genotypes[0]["NA00001"] == "0/1"
    assert genotypes[1]["SAMPLE2"] == "./."

    memberships = read_tsv(out_dir / "cluster_memberships.chr1.tsv")
    assert len(memberships) == 4
    assert {(row["sample_id"], row["haplotype"], row["cluster_id"]) for row in memberships} == {
        ("NA00001", "0", "cluster_a"),
        ("NA00001", "1", "cluster_b"),
        ("SAMPLE2", "0", "cluster_a"),
        ("SAMPLE2", "1", "cluster_b"),
    }

    block_summary = read_tsv(out_dir / "sv_block_summary.chr1.tsv")
    assert len(block_summary) == 2
    assert {(row["sv_id"], row["haploblock_id"]) for row in block_summary} == {
        ("sv1", "chr1_0_1000"),
        ("sv2", "chr1_0_1000"),
    }


def test_legacy_cleanup_keeps_published_contracts(tmp_path: Path) -> None:
    membership_path = tmp_path / "cluster_memberships.chr6.tsv"
    membership_path.write_text("published\n")
    legacy_path = tmp_path / "cluster_evidence.tsv"
    legacy_path.write_text("obsolete\n")

    match_svs_to_clusters.remove_legacy_outputs(tmp_path)

    assert membership_path.read_text() == "published\n"
    assert not legacy_path.exists()
