"""Data contracts for Problem A (herb drying): inputs and result workbooks."""

from __future__ import annotations

from forge.contracts import Column, FrameContract
from forge.xlsx import SheetContract, WorkbookContract

INPUT_CONTRACTS: dict[str, FrameContract] = {
    "附件1__Sheet1": FrameContract(
        name="A.附件1.烘房温湿度",
        columns=(
            Column("时间", "int", min=0, unique=True, monotonic="increasing", step=60),
            Column("温度", "float", min=0, max=120),
            Column("水分浓度", "float", min=0, max=1),
        ),
        min_rows=200,
    ),
    "附件2__Sheet1": FrameContract(
        name="A.附件2.药材半径",
        columns=(
            Column("时间", "int", min=0, unique=True, monotonic="increasing"),
            Column("半径", "float", min=0.1, max=2.0, monotonic="non_increasing"),
        ),
        min_rows=100,
    ),
}

_RADIAL = 22  # A1 label + 21 distances 0.0 … 2.0 (step 0.1 cm)
_HEADER = [0.0, 0.5, 1.0, 1.5, 2.0]  # first-row distances checked at columns B, G, L, Q, V (A1 is the label)
WHOLE_PROCESS_FILE = "result2_全过程.xlsx"  # supplementary deliverable, whole process at 60 s (MDR-0005)


def _sheet(name: str, step: float, **kw: object) -> SheetContract:
    """Result sheet contract: header positions, time-column step and the numeric region."""
    return SheetContract(
        name,
        header_len=kw.pop("header_len", _RADIAL),  # type: ignore[arg-type]
        numeric_from_col=0,
        first_col_step=step,
        header_values=tuple(_HEADER),
        header_value_columns=(1, 6, 11, 16, 21),
        **kw,  # type: ignore[arg-type]
    )


RESULT_CONTRACTS: list[WorkbookContract] = [
    WorkbookContract(
        "result1.xlsx",
        "result1.xlsx",
        (
            _sheet("温度", 1.0, min_rows=1800, max_rows=1801),
            _sheet("水分浓度", 1.0, min_rows=1800, max_rows=1801),
        ),
    ),
    WorkbookContract(
        "result2.xlsx",
        "result2.xlsx",
        (
            _sheet("温度", 1.0, min_rows=10800, max_rows=10801),
            _sheet("水分浓度", 1.0, min_rows=10800, max_rows=10801),
        ),
    ),
    WorkbookContract("result3.xlsx", "result3.xlsx", (_sheet("Sheet1", 60.0, min_rows=60),)),
    WorkbookContract(
        "result4.xlsx",
        "result4.xlsx",
        (_sheet("Sheet1", 60.0, min_rows=60, header_len=_RADIAL + 1, allow_blank=True),),
    ),
    WorkbookContract(
        WHOLE_PROCESS_FILE,
        "result2.xlsx",
        (_sheet("温度", 60.0, min_rows=60), _sheet("水分浓度", 60.0, min_rows=60)),
    ),
]
