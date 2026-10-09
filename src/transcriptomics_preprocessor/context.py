"""Objects shared by all stage implementations."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any

from .config import PipelineConfig
from .layout import OutputLayout

if TYPE_CHECKING:
    from .console import Console


@dataclass
class StageContext:
    """Everything a stage needs besides the sample itself."""

    layout: OutputLayout
    config: PipelineConfig
    tmp_root: Path
    console: Console
    tool_versions: dict[str, str | None]
    threads: int


@dataclass
class StageResult:
    """What a stage produced.

    ``outputs`` are files recorded in the run journal (and checked on resume);
    ``metrics`` are measured values; ``notes`` are human-readable statements
    shown after the stage completes; ``warnings`` require attention.
    """

    outputs: dict[str, Path] = field(default_factory=dict)
    metrics: dict[str, Any] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
