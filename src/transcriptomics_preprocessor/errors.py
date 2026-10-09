"""Exception hierarchy and process exit codes.

Every error raised deliberately by the application derives from
:class:`PreprocessorError`, carries a stable exit code and, where possible, a
``hint`` telling the user what to do next.
"""

from __future__ import annotations


class ExitCode:
    """Process exit codes used consistently by the command-line interface."""

    OK = 0
    #: Unexpected internal error (a bug); details are in the log.
    INTERNAL = 1
    #: Invalid command-line arguments, configuration or metadata.
    USAGE = 2
    #: A required external tool or Python dependency is missing or broken.
    DEPENDENCY = 3
    #: The run finished but at least one sample or report stage failed.
    PARTIAL = 4
    #: No sample completed successfully, or a run-level failure occurred.
    FAILED = 5
    #: Interrupted by the user (SIGINT) or terminated (SIGTERM).
    INTERRUPTED = 130


class PreprocessorError(Exception):
    """Base class for expected, user-facing errors."""

    exit_code: int = ExitCode.FAILED

    def __init__(self, message: str, hint: str | None = None) -> None:
        super().__init__(message)
        self.message = message
        self.hint = hint

    def __str__(self) -> str:
        return self.message


class InputValidationError(PreprocessorError):
    """Invalid user input: accessions, metadata, paths or parameters."""

    exit_code = ExitCode.USAGE

    def __init__(
        self, message: str, hint: str | None = None, problems: list[str] | None = None
    ) -> None:
        super().__init__(message, hint)
        self.problems: list[str] = list(problems or [])


class ConfigError(InputValidationError):
    """Invalid configuration file or configuration value."""


class DependencyError(PreprocessorError):
    """A required external program is missing or unusable."""

    exit_code = ExitCode.DEPENDENCY


class OutputDirectoryError(PreprocessorError):
    """The output directory cannot be used safely as requested."""

    exit_code = ExitCode.USAGE


class StateError(PreprocessorError):
    """The run journal (checkpoint file) is unreadable or inconsistent."""

    exit_code = ExitCode.FAILED


class StageError(PreprocessorError):
    """A processing stage failed for one sample or one aggregate step."""

    exit_code = ExitCode.FAILED


class ToolExecutionError(StageError):
    """An external program exited unsuccessfully or timed out."""

    def __init__(
        self,
        message: str,
        hint: str | None = None,
        returncode: int | None = None,
        log_path: str | None = None,
        output_tail: str = "",
    ) -> None:
        super().__init__(message, hint)
        self.returncode = returncode
        self.log_path = log_path
        self.output_tail = output_tail


class IntegrityError(StageError):
    """An expected output is missing, empty, truncated or malformed."""


class DiskSpaceError(StageError):
    """Not enough free disk space for the next step."""


class InterruptedRun(PreprocessorError):
    """The run was interrupted by a signal."""

    exit_code = ExitCode.INTERRUPTED
