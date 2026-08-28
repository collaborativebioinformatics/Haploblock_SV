"""Stage 9: summarize Stages 5-8 and write a biological interpretation report.

Statistics and figures are generated locally from the stage tables. When an
OpenAI API key is available, the compact ``report_facts.json`` is also sent to
the Responses API for a plain-language interpretation. Raw genotype and
association tables are never sent to the model.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import matplotlib
import mistune

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


AGENT_INSTRUCTIONS = """You are a careful population-genetics reporting assistant.
Write a concise biological interpretation of the supplied pipeline facts in Markdown.
Use the sections: Main findings, Biological interpretation, Limitations, and Next steps.
Every number must come directly from the supplied JSON. Distinguish statistical results
from biological hypotheses. Gene overlap is annotation, not evidence of functional
impact. PCA is quality control, not proof that haploblocks explain population structure.
Do not add literature claims or outside facts. Explicitly mention when no associations
pass the reported q-value threshold or when reported permutation counts limit resolution.
"""


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config", type=Path, action="append", required=True,
        help=(
            "Stage config containing registered output paths. Repeat to combine "
            "branches, normally once for Stage 5 and once for Stage 8."
        ),
    )
    parser.add_argument("--out-dir", type=Path, default=Path("stage9_output"))
    parser.add_argument(
        "--agent", choices=("auto", "required", "off"), default="auto",
        help=(
            "auto uses the agent when OPENAI_API_KEY is set; required fails without "
            "a key; off writes only the deterministic report"
        ),
    )
    parser.add_argument("--model", default="gpt-5.6")
    parser.add_argument("--top-candidates", type=int, default=10)
    return parser.parse_args(argv)


def resolved_paths(config_paths: list[Path]) -> dict[str, Path]:
    """Merge scalar paths from carried-forward stage configs."""
    paths: dict[str, Path] = {}
    for config_path in config_paths:
        config = yaml.safe_load(config_path.read_text())
        for name, value in config["paths"].items():
            if isinstance(value, str):
                path = Path(value)
                paths[name] = path if path.is_absolute() else config_path.parent / path
    return paths


def clean_value(value: object) -> object:
    if pd.isna(value):
        return None
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return round(float(value), 6)
    return str(value)


def records(table: pd.DataFrame, columns: list[str], limit: int) -> list[dict[str, object]]:
    present = [column for column in columns if column in table]
    return [
        {column: clean_value(value) for column, value in row.items()}
        for row in table[present].head(limit).to_dict("records")
    ]


def value_counts(table: pd.DataFrame, column: str) -> dict[str, int]:
    if column not in table:
        return {}
    return {
        str(name): int(count)
        for name, count in table[column].fillna("not_reported").value_counts().items()
    }


def build_facts(
    enrichment: pd.DataFrame,
    associations: pd.DataFrame,
    annotations: pd.DataFrame,
    representation: pd.DataFrame | None,
    information: pd.DataFrame | None,
    pca_variance: pd.DataFrame | None,
    top_candidates: int,
) -> dict[str, object]:
    enrichment = enrichment.copy()
    enrichment["fold_enrichment"] = np.divide(
        enrichment["observed_count"], enrichment["expected_count"],
        out=np.full(len(enrichment), np.nan),
        where=enrichment["expected_count"].to_numpy() > 0,
    )
    stage5_top = enrichment.sort_values(
        ["q_value", "fold_enrichment"], ascending=[True, False]
    )

    associations = associations.copy()
    associations["absolute_population_adjusted_r"] = associations[
        "population_adjusted_r"
    ].abs()
    stage6_top = associations.sort_values(
        ["q_value", "absolute_population_adjusted_r"], ascending=[True, False]
    )
    significant = associations["q_value"] < 0.05

    annotations = annotations.copy()
    annotations["absolute_population_adjusted_r"] = annotations[
        "population_adjusted_r"
    ].abs()
    annotations["association_significant"] = annotations["q_value"].lt(0.05)
    if "call_quality" in annotations:
        annotations["pass_precise"] = annotations["call_quality"].eq("pass_precise")
    else:
        annotations["pass_precise"] = annotations["filter"].eq("PASS") & ~annotations[
            "imprecise"
        ].astype(bool)
    if "consequence_priority" not in annotations:
        annotations["consequence_priority"] = annotations.get("candidate_score", 0)
    supported_patterns = {
        "hash_tag_candidate", "multi_cluster_sv_candidate", "hash_subdivision_candidate",
        "multi_cluster_and_subdivision_candidate",
    }
    annotations["representation_supported"] = (
        annotations["representation_pattern"].isin(supported_patterns)
        if "representation_pattern" in annotations else False
    )
    stage8_top = annotations.sort_values(
        ["association_significant", "pass_precise", "consequence_priority",
         "representation_supported", "absolute_population_adjusted_r"],
        ascending=[False, False, False, False, False],
    )

    stage7: dict[str, object]
    if representation is not None:
        stage7 = {
            "n_sv_haploblock_pairs": int(len(representation)),
            "representation_pattern_counts": value_counts(
                representation, "representation_pattern"
            ),
            "population_context_pattern_counts": value_counts(
                representation, "population_context_pattern"
            ),
        }
    elif information is not None:
        stage7 = {
            "n_sv_haploblock_pairs": int(len(information)),
            "median_normalized_information_gain": clean_value(
                information["normalized_information_gain"].median()
            ),
            "pairs_with_information_gain_at_least_0_5": int(
                information["normalized_information_gain"].ge(0.5).sum()
            ),
            "median_mixed_diplotype_fraction": clean_value(
                information["mixed_diplotype_fraction"].median()
            ),
            "note": "Legacy Stage 7 output; direct representation categories were unavailable.",
        }
    else:
        stage7 = {"note": "Stage 7 summary was not registered in the supplied configs."}

    if pca_variance is not None:
        stage7["pca"] = {
            "n_variants": int(pca_variance["n_variants"].max()),
            "explained_variance": records(
                pca_variance, ["component", "explained_variance_ratio"], len(pca_variance)
            ),
            "interpretation": "QC only",
        }

    min_p = associations["p_value"].min()
    permutation_counts = (
        value_counts(associations, "permutations_used")
        if "permutations_used" in associations else {}
    )
    return {
        "stage5_sv_type_enrichment": {
            "n_block_type_cells": int(len(enrichment)),
            "n_significant_q_lt_0_05": int(enrichment["q_value"].lt(0.05).sum()),
            "n_flagged_by_stage": int(enrichment["flagged"].sum()),
            "sv_types": sorted(enrichment["sv_type"].astype(str).unique()),
            "top_cells": records(
                stage5_top,
                ["haploblock_id", "sv_type", "observed_count", "expected_count",
                 "fold_enrichment", "p_value", "q_value", "flagged"],
                top_candidates,
            ),
            "model_limitation": "Expected counts are adjusted for haploblock length only.",
        },
        "stage6_cluster_association": {
            "n_sv_haploblock_pairs": int(len(associations)),
            "n_significant_q_lt_0_05": int(significant.sum()),
            "minimum_p_value": clean_value(min_p),
            "permutations_used_counts": permutation_counts,
            "association_pattern_counts": value_counts(associations, "association_pattern"),
            "effect_size_counts": {
                f"absolute_r_at_least_{threshold}": int(
                    associations["absolute_population_adjusted_r"].ge(threshold).sum()
                )
                for threshold in (0.3, 0.5, 0.7, 0.9)
            },
            "top_associations": records(
                stage6_top,
                ["sv_record_id", "sv_id", "chrom", "start", "end", "sv_type",
                 "haploblock_id", "best_cluster_id", "population_adjusted_r", "p_value",
                 "q_value", "carrier_rate_with_cluster", "carrier_rate_without_cluster",
                 "association_pattern"],
                top_candidates,
            ),
        },
        "stage7_hash_representation_and_qc": stage7,
        "stage8_candidate_annotation": {
            "n_annotated_sv_haploblock_pairs": int(len(annotations)),
            "consequence_counts": value_counts(annotations, "consequence"),
            "sv_type_counts": value_counts(annotations, "sv_type"),
            "top_candidates": records(
                stage8_top,
                ["sv_record_id", "sv_id", "chrom", "start", "end", "sv_type",
                 "haploblock_id", "genes", "consequence", "overlap_basis", "filter",
                 "imprecise", "population_adjusted_r", "q_value", "association_pattern",
                 "representation_pattern", "sv_class", "specific_to_population"],
                top_candidates,
            ),
            "candidate_selection": (
                "Ordered by association significance, precise PASS status, consequence priority, "
                "reported representation support, then absolute population-adjusted correlation."
            ),
            "interpretation_limit": (
                "Coordinate-based consequence labels nominate candidates but do not establish "
                "molecular function."
            ),
        },
    }


def save_bar_plot(counts: dict[str, int], title: str, path: Path) -> bool:
    if not counts:
        return False
    items = sorted(counts.items(), key=lambda item: item[1], reverse=True)[:12][::-1]
    labels, values = zip(*items)
    fig, ax = plt.subplots(figsize=(8, max(3.5, 0.38 * len(items))))
    ax.barh(labels, values, color="#4C78A8")
    ax.set_title(title)
    ax.set_xlabel("Count")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def save_stage5_plot(top_cells: list[dict[str, object]], path: Path) -> bool:
    rows = [
        row for row in top_cells
        if row["fold_enrichment"] is not None and float(row["fold_enrichment"]) > 0
    ]
    if not rows:
        return False
    rows = rows[::-1]
    labels = [f"{row['haploblock_id']} · {row['sv_type']}" for row in rows]
    values = [np.log2(float(row["fold_enrichment"])) for row in rows]
    colors = ["#E45756" if value >= 0 else "#4C78A8" for value in values]
    fig, ax = plt.subplots(figsize=(9, max(4, 0.42 * len(rows))))
    ax.barh(labels, values, color=colors)
    ax.axvline(0, color="#666666", linewidth=0.8)
    ax.set_title("Stage 5 top block/type enrichment results")
    ax.set_xlabel("log2(observed / expected)")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return True


def create_figures(
    facts: dict[str, object],
    information: pd.DataFrame | None,
    figure_dir: Path,
) -> list[str]:
    figure_dir.mkdir(parents=True, exist_ok=True)
    figures: list[str] = []
    stage5_filename = "stage5_type_enrichment.png"
    if save_stage5_plot(
        facts["stage5_sv_type_enrichment"]["top_cells"], figure_dir / stage5_filename
    ):
        figures.append(f"figures/{stage5_filename}")
    plot_specs = [
        (
            facts["stage6_cluster_association"]["association_pattern_counts"],
            "Stage 6 association patterns", "stage6_association_patterns.png",
        ),
        (
            facts["stage8_candidate_annotation"]["consequence_counts"],
            "Stage 8 consequence annotations", "stage8_consequences.png",
        ),
    ]
    stage7 = facts["stage7_hash_representation_and_qc"]
    if "representation_pattern_counts" in stage7:
        plot_specs.append((
            stage7["representation_pattern_counts"],
            "Stage 7 representation patterns", "stage7_representation_patterns.png",
        ))
    for counts, title, filename in plot_specs:
        if save_bar_plot(counts, title, figure_dir / filename):
            figures.append(f"figures/{filename}")

    if information is not None and "representation_pattern_counts" not in stage7:
        values = information["normalized_information_gain"].dropna()
        if not values.empty:
            fig, ax = plt.subplots(figsize=(7, 4))
            ax.hist(values, bins=30, color="#4C78A8", edgecolor="white")
            ax.set(title="Stage 7 normalized information gain", xlabel="Normalized information gain",
                   ylabel="SV-haploblock pairs")
            ax.spines[["top", "right"]].set_visible(False)
            fig.tight_layout()
            filename = "stage7_information_gain.png"
            fig.savefig(figure_dir / filename, dpi=160)
            plt.close(fig)
            figures.append(f"figures/{filename}")
    return figures


def request_agent_interpretation(
    facts: dict[str, object], model: str, client: object | None = None
) -> tuple[str, str | None]:
    if client is None:
        from openai import OpenAI

        client = OpenAI()
    response = client.responses.create(
        model=model,
        instructions=AGENT_INSTRUCTIONS,
        input=json.dumps(facts, indent=2, allow_nan=False),
    )
    return response.output_text.strip(), getattr(response, "id", None)


def factual_markdown(facts: dict[str, object], figures: list[str]) -> str:
    stage5 = facts["stage5_sv_type_enrichment"]
    stage6 = facts["stage6_cluster_association"]
    stage7 = facts["stage7_hash_representation_and_qc"]
    stage8 = facts["stage8_candidate_annotation"]
    lines = [
        "# Stage 5-8 biological report", "", "## Deterministic summary", "",
        f"- Stage 5 tested {stage5['n_block_type_cells']:,} block/type cells; "
        f"{stage5['n_significant_q_lt_0_05']:,} had q < 0.05.",
        f"- Stage 6 summarized {stage6['n_sv_haploblock_pairs']:,} SV-haploblock pairs; "
        f"{stage6['n_significant_q_lt_0_05']:,} passed q < 0.05.",
        f"- Stage 7 summarized {stage7.get('n_sv_haploblock_pairs', 0):,} "
        "SV-haploblock pairs where a registered table was available.",
        f"- Stage 8 annotated {stage8['n_annotated_sv_haploblock_pairs']:,} "
        "SV-haploblock pairs.", "",
        "The exact values used for interpretation are saved in `report_facts.json`.", "",
        "## Figures", "",
    ]
    for figure in figures:
        lines.extend([f"![{Path(figure).stem.replace('_', ' ')}]({figure})", ""])
    lines.extend([
        "## Required interpretation cautions", "",
        f"- {stage5['model_limitation']}",
        "- Stage 6 effect sizes are candidates until supported by calibrated permutation "
        "tests, technical-batch controls, and validation across populations.",
        "- Stage 7 PCA is a quality-control view, not evidence that haploblocks explain ancestry.",
        f"- {stage8['interpretation_limit']}", "",
    ])
    return "\n".join(lines)


def write_html(markdown_text: str, path: Path) -> None:
    body = mistune.html(markdown_text)
    path.write_text(f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>Stage 5-8 biological report</title>
<style>body{{max-width:900px;margin:40px auto;padding:0 20px;font:16px/1.55 system-ui;color:#202124}}
h1,h2,h3{{line-height:1.2}} img{{max-width:100%;height:auto}} code{{background:#f2f2f2;padding:2px 4px}}
li{{margin:.35rem 0}}</style></head><body>{body}</body></html>\n""")


def run_report(
    config_paths: list[Path],
    out_dir: Path,
    agent_mode: str = "auto",
    model: str = "gpt-5.6",
    top_candidates: int = 10,
) -> None:
    paths = resolved_paths(config_paths)
    enrichment = pd.read_csv(paths["sv_type_enrichment"], sep="\t")
    associations = pd.read_csv(paths["sv_cluster_summary"], sep="\t")
    annotations = pd.read_csv(paths["annotated_sv_candidates"], sep="\t")
    representation = (
        pd.read_csv(paths["sv_hash_representation"], sep="\t")
        if "sv_hash_representation" in paths else None
    )
    information = (
        pd.read_csv(paths["sv_haploblock_information"], sep="\t")
        if "sv_haploblock_information" in paths else None
    )
    pca_variance = (
        pd.read_csv(paths["sv_pca_variance"], sep="\t")
        if "sv_pca_variance" in paths else None
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    facts = build_facts(
        enrichment, associations, annotations, representation, information,
        pca_variance, top_candidates,
    )
    (out_dir / "report_facts.json").write_text(
        json.dumps(facts, indent=2, allow_nan=False) + "\n"
    )
    figures = create_figures(facts, information, out_dir / "figures")
    report = factual_markdown(facts, figures)

    api_key_available = bool(os.environ.get("OPENAI_API_KEY"))
    use_agent = agent_mode == "required" or (
        agent_mode == "auto" and api_key_available
    )
    metadata: dict[str, object] = {
        "agent_mode": agent_mode, "model": model, "agent_used": use_agent,
        "prompt_version": 1,
        "response_id": None,
    }
    agent_error = (
        "RuntimeError: OPENAI_API_KEY is required when --agent required is used"
        if agent_mode == "required" and not api_key_available else None
    )
    if agent_error is None and use_agent:
        try:
            interpretation, response_id = request_agent_interpretation(facts, model)
            report += "\n## Agent interpretation\n\n" + interpretation + "\n"
            metadata["response_id"] = response_id
        except Exception as error:
            agent_error = f"{type(error).__name__}: {error}"
    if agent_error is not None:
        metadata["agent_used"] = False
        metadata["error"] = agent_error
        report += (
            "\n## Agent interpretation\n\nNot generated because the agent request failed. "
            "The deterministic report is complete; see `agent_metadata.json` for details.\n"
        )
    elif not use_agent:
        report += (
            "\n## Agent interpretation\n\nNot generated. Set `OPENAI_API_KEY` and use "
            "`--agent auto` or `--agent required` to add it.\n"
        )
    (out_dir / "agent_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    (out_dir / "report.md").write_text(report)
    write_html(report, out_dir / "report.html")
    if agent_error is not None and agent_mode == "required":
        raise RuntimeError(f"Agent interpretation failed: {agent_error}")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    run_report(args.config, args.out_dir, args.agent, args.model, args.top_candidates)


if __name__ == "__main__":
    main(sys.argv[1:])
