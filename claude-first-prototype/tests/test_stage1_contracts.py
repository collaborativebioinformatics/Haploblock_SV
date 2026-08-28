"""Contract tests for Stage 1's downstream tables."""

import csv
import gzip
import json
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
            "--gtf", str(tmp_path / "genes.gtf"),
            "--cluster-root", str(tmp_path / "clusters"),
            "--chroms", "chr1",
            "--threads", "1",
            "--out-dir", str(out_dir),
        ]
    )

    config = yaml.safe_load((out_dir / "config.yaml").read_text())
    assert {
        "vcf",
        "samples",
        "sample_metadata",
        "sv_genotypes",
        "sv_to_clusters",
        "haploblocks",
        "cluster_memberships",
        "sv_block_summary",
        "debug_and_qc",
    }.issubset(config["paths"])
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
    assert [row["sv_record_id"] for row in genotypes] == ["chr1_record_1", "chr1_record_2"]
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


def test_stage1_downloads_default_ont_metadata(tmp_path: Path, monkeypatch) -> None:
    def write_samples(args) -> None:
        args.out_dir.mkdir(parents=True, exist_ok=True)
        (args.out_dir / "samples.tsv").write_text(
            "sample_id\toriginal_sample_id\nNA00001\tGM00001\nSAMPLE2\tSAMPLE2\n"
        )

    downloaded_metadata = (
        "NHGRI_ID\tSex\tSubPopulation\tSuperPopulation\tONT_library\tONT_pore\n"
        "GM00001\tXY\tPOP1\tSUPER1\tLSK110\tR9\n"
        "SAMPLE2\tXX\tPOP2\tSUPER2\tLSK114\tR10\n"
    ).encode()
    monkeypatch.setattr(stage1_cluster_aware, "run_cluster_aware", write_samples)
    monkeypatch.setattr(
        stage1_cluster_aware,
        "request_with_retries",
        lambda url, retries: downloaded_metadata,
    )

    out_dir = tmp_path / "stage1_output"
    stage1_cluster_aware.main(
        ["--vcf", str(tmp_path / "cohort.vcf.gz"), "--chroms", "chr1", "--out-dir", str(out_dir)]
    )

    assert read_tsv(out_dir / "sample_metadata.tsv") == [
        {
            "sample_id": "NA00001",
            "Sex": "XY",
            "population": "POP1",
            "superpopulation": "SUPER1",
            "ONT_library": "LSK110",
            "ONT_pore": "R9",
        },
        {
            "sample_id": "SAMPLE2",
            "Sex": "XX",
            "population": "POP2",
            "superpopulation": "SUPER2",
            "ONT_library": "LSK114",
            "ONT_pore": "R10",
        },
    ]
    config = yaml.safe_load((out_dir / "config.yaml").read_text())
    assert config["data_sources"]["sample_metadata"] == (
        stage1_cluster_aware.DEFAULT_SAMPLE_METADATA_URL
    )
    assert config["paths"]["sample_metadata"] == str(
        (out_dir / "sample_metadata.tsv").resolve()
    )
    assert config["data_sources"]["gtf"] == stage1_cluster_aware.DEFAULT_GTF_URL
    assert config["paths"]["gtf"] == str(
        (out_dir / "Homo_sapiens.GRCh38.115.gtf.gz").resolve()
    )


def test_legacy_cleanup_keeps_published_contracts(tmp_path: Path) -> None:
    membership_path = tmp_path / "cluster_memberships.chr6.tsv"
    membership_path.write_text("published\n")
    legacy_path = tmp_path / "cluster_evidence.tsv"
    legacy_path.write_text("obsolete\n")

    match_svs_to_clusters.remove_legacy_outputs(tmp_path)

    assert membership_path.read_text() == "published\n"
    assert not legacy_path.exists()


def test_downloaded_cluster_files_are_reused(tmp_path: Path, monkeypatch) -> None:
    names = ["chr1_0-100_cluster.tsv", "chr1_100-200_cluster.tsv"]
    downloads = []
    monkeypatch.setattr(
        match_svs_to_clusters,
        "discover_cluster_filenames",
        lambda base_url, chrom, retries: names,
    )

    def download(url: str, retries: int) -> bytes:
        downloads.append(url)
        return b"cluster_a\tsample_chr1_region_0-100_hap0\n"

    monkeypatch.setattr(match_svs_to_clusters, "request_with_retries", download)
    cache_dir = match_svs_to_clusters.cluster_cache_directory(
        tmp_path, "https://example.test/source", "chr1"
    )
    assert cache_dir == match_svs_to_clusters.cluster_cache_directory(
        tmp_path, "https://example.test/source/", "chr1"
    )
    assert cache_dir != match_svs_to_clusters.cluster_cache_directory(
        tmp_path, "https://other.example.test/source", "chr1"
    )
    first_paths, first_qc = match_svs_to_clusters.prepare_cluster_files(
        "chr1", "https://example.test", cache_dir, None, 2, 1
    )
    second_paths, second_qc = match_svs_to_clusters.prepare_cluster_files(
        "chr1", "https://example.test", cache_dir, None, 2, 1
    )

    assert first_paths == second_paths
    assert len(downloads) == 2
    assert first_qc["downloaded"] == 2
    assert second_qc["downloaded"] == 0
    assert second_qc["reused"] == 2


def test_parallel_block_evaluation_preserves_outputs(tmp_path: Path) -> None:
    vcf_path = tmp_path / "cohort.vcf.gz"
    write_test_vcf(vcf_path)
    cluster_dir = tmp_path / "clusters" / "chr1"
    cluster_dir.mkdir(parents=True)
    membership = (
        "cluster_a\tGM00001_chr1_region_{start}-{end}_hap0\n"
        "cluster_b\tGM00001_chr1_region_{start}-{end}_hap1\n"
        "cluster_a\tSAMPLE2_chr1_region_{start}-{end}_hap0\n"
        "cluster_b\tSAMPLE2_chr1_region_{start}-{end}_hap1\n"
    )
    for start, end in [(0, 1000), (50, 500)]:
        (cluster_dir / f"chr1_{start}-{end}_cluster.tsv").write_text(
            membership.format(start=start, end=end)
        )

    output_dirs = []
    for threads in [1, 2]:
        out_dir = tmp_path / f"output_{threads}"
        args = match_svs_to_clusters.parse_args(
            [
                "--vcf", str(vcf_path),
                "--chroms", "chr1",
                "--cluster-root", str(tmp_path / "clusters"),
                "--out-dir", str(out_dir),
                "--threads", str(threads),
            ]
        )
        match_svs_to_clusters.run(args)
        output_dirs.append(out_dir)

    for name in [
        "sv_genotypes.chr1.tsv",
        "haploblocks.chr1.tsv",
        "cluster_memberships.chr1.tsv",
        "sv_block_summary.chr1.tsv",
        "sv_to_clusters.chr1.tsv",
    ]:
        assert (output_dirs[0] / name).read_text() == (output_dirs[1] / name).read_text()

    qc_fields = [
        "sv_block_pairs",
        "converged_models",
        "nonconverged_models",
        "cluster_evaluations",
        "associated_cluster_rows",
        "cluster_probability_bins",
        "cluster_call_rate_bins",
        "heterozygote_assignment_status_counts",
        "heterozygote_assignment_confidence_bins",
        "association_class_counts",
    ]
    qcs = [
        json.loads((out_dir / "debug_and_qc" / "method_qc.chr1.json").read_text())
        for out_dir in output_dirs
    ]
    assert {field: qcs[0][field] for field in qc_fields} == {
        field: qcs[1][field] for field in qc_fields
    }


def test_single_worker_processes_every_requested_chromosome(tmp_path: Path, monkeypatch) -> None:
    out_dir = tmp_path / "output"
    (out_dir / "debug_and_qc").mkdir(parents=True)
    vcf_path = tmp_path / "cohort.vcf.gz"
    processed = []
    monkeypatch.setattr(match_svs_to_clusters, "chromosome_is_complete", lambda *args: False)
    monkeypatch.setattr(
        match_svs_to_clusters,
        "reusable_sv_tables",
        lambda vcf, chroms, output, max_length: {
            chrom: tmp_path / f"{chrom}.tsv" for chrom in chroms
        },
    )
    monkeypatch.setattr(
        match_svs_to_clusters,
        "process_chromosome",
        lambda chrom, *args: processed.append(chrom),
    )
    args = match_svs_to_clusters.parse_args(
        [
            "--vcf", str(vcf_path),
            "--chroms", "chr1,chr2",
            "--out-dir", str(out_dir),
            "--threads", "1",
        ]
    )

    match_svs_to_clusters.run(args)

    assert processed == ["chr1", "chr2"]


def test_completed_chromosome_must_match_cluster_source(tmp_path: Path) -> None:
    out_dir = tmp_path / "output"
    qc_dir = out_dir / "debug_and_qc"
    qc_dir.mkdir(parents=True)
    vcf_path = tmp_path / "cohort.vcf.gz"
    vcf_path.touch()
    for path in [
        out_dir / "samples.tsv",
        out_dir / "sv_genotypes.chr1.tsv",
        out_dir / "sv_to_clusters.chr1.tsv",
        out_dir / "haploblocks.chr1.tsv",
        out_dir / "cluster_memberships.chr1.tsv",
        out_dir / "sv_block_summary.chr1.tsv",
    ]:
        path.touch()
    (qc_dir / "vcf_qc.chr1.json").write_text(
        json.dumps({"vcf": str(vcf_path.resolve()), "sv_id_max_length": 80})
    )
    (qc_dir / "method_qc.chr1.json").write_text(
        json.dumps(
            {
                "cluster_download": {
                    "source": "download",
                    "base_url": "https://source-a.example.test",
                },
                "association_threshold": 0.75,
                "posterior_assignment_threshold": 0.75,
                "max_iterations": 25,
                "tolerance": 1e-5,
            }
        )
    )
    common = ("chr1", vcf_path, out_dir)
    settings = (80, 0.75, 0.75, 25, 1e-5)

    assert match_svs_to_clusters.chromosome_is_complete(
        *common, "https://source-a.example.test/", None, *settings
    )
    assert not match_svs_to_clusters.chromosome_is_complete(
        *common, "https://source-b.example.test", None, *settings
    )
