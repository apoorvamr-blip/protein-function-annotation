"""Runtime configuration.

Configuration is assembled from (in increasing priority) built-in defaults,
environment variables (optionally loaded from a ``.env`` file) and explicit
overrides supplied by the CLI. No credentials or personal data are hard-coded.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

from protein_annotator import __version__

try:  # python-dotenv is a declared dependency, but keep config importable without it.
    from dotenv import find_dotenv, load_dotenv
except ImportError:  # pragma: no cover
    find_dotenv = None  # type: ignore[assignment]
    load_dotenv = None  # type: ignore[assignment]

#: NCBI asks clients not to poll a single RID more than once per minute.
NCBI_MIN_POLL_INTERVAL = 60.0


class ConfigurationError(ValueError):
    """Raised when a configuration value is invalid."""


def _env_str(name: str) -> str | None:
    value = os.environ.get(name)
    if value is None:
        return None
    value = value.strip()
    return value or None


def _env_float(name: str, default: float) -> float:
    value = _env_str(name)
    if value is None:
        return default
    try:
        return float(value)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be a number, got {value!r}") from exc


def _env_int(name: str, default: int) -> int:
    value = _env_str(name)
    if value is None:
        return default
    try:
        return int(value)
    except ValueError as exc:
        raise ConfigurationError(f"{name} must be an integer, got {value!r}") from exc


def _env_bool(name: str, default: bool) -> bool:
    value = _env_str(name)
    if value is None:
        return default
    lowered = value.lower()
    if lowered in {"1", "true", "yes", "on"}:
        return True
    if lowered in {"0", "false", "no", "off"}:
        return False
    raise ConfigurationError(f"{name} must be a boolean, got {value!r}")


def default_cache_dir() -> Path:
    """Return the default on-disk cache directory (XDG compliant)."""
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "protein-annotator"


@dataclass(frozen=True)
class Endpoints:
    """Base URLs of the remote services."""

    blast_url: str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
    interpro_url: str = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
    uniprot_url: str = "https://rest.uniprot.org"


@dataclass(frozen=True)
class BlastSettings:
    """NCBI BLASTP (URL API) settings."""

    database: str = "swissprot"
    evalue: float = 1e-5
    max_hits: int = 25
    timeout: float = 1800.0
    poll_interval: float = NCBI_MIN_POLL_INTERVAL


@dataclass(frozen=True)
class InterProSettings:
    """EMBL-EBI InterProScan 5 job-dispatcher settings."""

    timeout: float = 1800.0
    poll_interval: float = 15.0
    applications: tuple[str, ...] = ()


@dataclass(frozen=True)
class UniProtSettings:
    """UniProt REST settings."""

    timeout: float = 60.0
    max_entries: int = 15
    map_non_uniprot_ids: bool = True
    mapping_timeout: float = 300.0
    mapping_poll_interval: float = 5.0


@dataclass(frozen=True)
class InferenceSettings:
    """Thresholds used by evidence integration (documented in docs/methodology.md)."""

    significance_evalue: float = 1e-5
    min_query_coverage: float = 0.5
    cluster_similarity: float = 0.5
    competing_similarity: float = 0.2
    competing_score_ratio: float = 0.5
    go_min_fraction: float = 0.5


@dataclass(frozen=True)
class AnnotatorConfig:
    """Complete pipeline configuration."""

    email: str | None = None
    tool_name: str = "protein-function-annotator"
    request_timeout: float = 60.0
    max_retries: int = 5
    backoff_factor: float = 2.0
    cache_dir: Path = field(default_factory=default_cache_dir)
    use_cache: bool = True
    cache_ttl_days: float = 30.0
    endpoints: Endpoints = field(default_factory=Endpoints)
    blast: BlastSettings = field(default_factory=BlastSettings)
    interpro: InterProSettings = field(default_factory=InterProSettings)
    uniprot: UniProtSettings = field(default_factory=UniProtSettings)
    inference: InferenceSettings = field(default_factory=InferenceSettings)

    @property
    def user_agent(self) -> str:
        """User-Agent header identifying this tool (and contact e-mail when configured)."""
        contact = f"; mailto:{self.email}" if self.email else ""
        return f"{self.tool_name}/{__version__} (python-requests{contact})"

    def with_overrides(self, **changes: Any) -> AnnotatorConfig:
        """Return a copy with top-level fields replaced."""
        return replace(self, **changes)

    def validate(self) -> None:
        """Check values that would otherwise cause confusing downstream errors."""
        if self.blast.max_hits < 1:
            raise ConfigurationError("max hits must be at least 1")
        if self.blast.evalue <= 0:
            raise ConfigurationError("E-value threshold must be positive")
        for label, value in (
            ("BLAST timeout", self.blast.timeout),
            ("InterPro timeout", self.interpro.timeout),
            ("UniProt timeout", self.uniprot.timeout),
            ("HTTP timeout", self.request_timeout),
        ):
            if value <= 0:
                raise ConfigurationError(f"{label} must be positive")
        if self.email is not None and "@" not in self.email:
            raise ConfigurationError(f"Invalid contact e-mail address: {self.email!r}")

    @classmethod
    def from_env(cls, env_file: str | Path | None = None) -> AnnotatorConfig:
        """Build configuration from environment variables (and an optional .env file)."""
        if load_dotenv is not None:
            if env_file is not None:
                load_dotenv(env_file, override=False)
            elif find_dotenv is not None:
                found = find_dotenv(usecwd=True)
                if found:
                    load_dotenv(found, override=False)

        email = _env_str("ANNOTATOR_EMAIL") or _env_str("NCBI_EMAIL") or _env_str("EBI_EMAIL")
        cache_dir_value = _env_str("ANNOTATOR_CACHE_DIR")
        applications = tuple(
            part.strip()
            for part in (_env_str("INTERPRO_APPLICATIONS") or "").split(",")
            if part.strip()
        )

        config = cls(
            email=email,
            request_timeout=_env_float("HTTP_TIMEOUT", 60.0),
            max_retries=_env_int("HTTP_MAX_RETRIES", 5),
            cache_dir=Path(cache_dir_value).expanduser() if cache_dir_value else default_cache_dir(),
            use_cache=_env_bool("ANNOTATOR_USE_CACHE", True),
            cache_ttl_days=_env_float("ANNOTATOR_CACHE_TTL_DAYS", 30.0),
            blast=BlastSettings(
                database=_env_str("BLAST_DATABASE") or "swissprot",
                evalue=_env_float("BLAST_EVALUE", 1e-5),
                max_hits=_env_int("BLAST_MAX_HITS", 25),
                timeout=_env_float("BLAST_TIMEOUT", 1800.0),
                poll_interval=max(
                    _env_float("BLAST_POLL_INTERVAL", NCBI_MIN_POLL_INTERVAL),
                    NCBI_MIN_POLL_INTERVAL,
                ),
            ),
            interpro=InterProSettings(
                timeout=_env_float("INTERPRO_TIMEOUT", 1800.0),
                poll_interval=max(_env_float("INTERPRO_POLL_INTERVAL", 15.0), 5.0),
                applications=applications,
            ),
            uniprot=UniProtSettings(
                timeout=_env_float("UNIPROT_TIMEOUT", 60.0),
                max_entries=_env_int("UNIPROT_MAX_ENTRIES", 15),
            ),
        )
        config.validate()
        return config
