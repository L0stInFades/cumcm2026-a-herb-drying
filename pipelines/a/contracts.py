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

RESULT_CONTRACTS: list[WorkbookContract] = [
    WorkbookContract("result1.xlsx", "result1.xlsx", (
        SheetContract("温度", min_rows=1800, max_rows=1801, header_len=_RADIAL, numeric_from_col=0),
        SheetContract("水分浓度", min_rows=1800, max_rows=1801, header_len=_RADIAL, numeric_from_col=0),
    )),
    WorkbookContract("result2.xlsx", "result2.xlsx", (
        SheetContract("温度", min_rows=10800, header_len=_RADIAL, numeric_from_col=0),
        SheetContract("水分浓度", min_rows=10800, header_len=_RADIAL, numeric_from_col=0),
    )),
    WorkbookContract("result3.xlsx", "result3.xlsx", (
        SheetContract("Sheet1", min_rows=60, header_len=_RADIAL, numeric_from_col=0),
    )),
    WorkbookContract("result4.xlsx", "result4.xlsx", (
        SheetContract("Sheet1", min_rows=60, numeric_from_col=0, allow_blank=True),
    )),
]
