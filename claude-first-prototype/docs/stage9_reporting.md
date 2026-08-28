# Stage 9: reporting

## Purpose

Stage 9 turns the selected Stage 5–8 outputs into a compact factual report. It is designed to present observed patterns, their limitations, and candidate loci without exposing full raw genotype and association tables to the reporting step.

## Inputs

The script receives one or more carried-forward `config.yaml` files. In the usual branch layout, supply the Stage 5 configuration for enrichment results and the Stage 8 configuration for association, representation, and annotation results. Registered scalar paths are merged from those configurations.

## Outputs

| File | Contents |
|---|---|
| `report_facts.json` | Compact counts, thresholds, categories, and top rows used to construct the report. |
| `report.md` | Deterministic Markdown findings and figures, with optional plain-language interpretation. |
| `report.html` | HTML rendering of the report. |
| `figures/` | Available enrichment, representation, and PCA figures. |
| `agent_metadata.json` | Status and metadata for optional agent interpretation. |

The report separates main findings, biological interpretation, limitations, and next steps. Report numbers derive from `report_facts.json`; that file is the auditable compact source for every summary statement.

## Optional agent interpretation

With `--agent auto` (the default), the script uses an OpenAI Responses API interpretation only when `OPENAI_API_KEY` is available. It receives `report_facts.json`, not raw genotype or full association data. `--agent off` writes only deterministic outputs; `--agent required` fails after writing deterministic artifacts if interpretation is unavailable.

## Running

```bash
python claude-first-prototype/pipeline/stage9_report.py \
  --config claude-first-prototype/stage5_output/config.yaml \
  --config claude-first-prototype/stage8_output/config.yaml \
  --out-dir claude-first-prototype/stage9_output \
  --agent off
```

The default report includes up to ten top candidates. The report prioritizes hypotheses for inspection; it does not replace validation, read-level review, or independent replication.
