"""Problem A stages: validate (here) plus the science and deliverable stages of the submodules.

Conventions (see docs/ENGINEERING_STANDARD.md):
  * every stage is ``@stage(name, deps=(...))`` and returns a JSON-serialisable metrics dict;
  * outputs go only to ``ctx.out(...)``; paper numbers via ``ctx.number("Key", value)``;
  * ``results`` writes result*.xlsx with ``forge.xlsx.write_result`` from the organisers' templates;
  * ``figures`` writes figures/*.pdf (+png preview) and ``tables`` writes tables/*.tex.

Stage map: q1, q2, q3, q4, convergence, verification, sensitivity (pipelines/a/science.py);
results, tables (pipelines/a/outputs.py); figures (pipelines/a/figures.py).
"""

from __future__ import annotations

from typing import Any

from forge.context import StageContext
from forge.runner import stage
from pipelines.a import extension, figures, outputs, science  # noqa: F401 - registers the stages on import
from pipelines.a.contracts import INPUT_CONTRACTS
from pipelines.common.validation import run_validation


@stage(
    "validate", deps=("ingest",), description="Validate the drying-chamber and radius series against their contracts"
)
def validate(ctx: StageContext) -> dict[str, Any]:
    return run_validation(ctx, INPUT_CONTRACTS)
