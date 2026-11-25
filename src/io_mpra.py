"""Compatibility imports for the legacy schema now named src.mpra_schema.

The context-preserving production reader is src.dataset, not this schema.
"""

from src.mpra_schema import load_mpra, tidy_mpra, validate_mpra

__all__ = ["load_mpra", "tidy_mpra", "validate_mpra"]
