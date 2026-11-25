"""Compatibility entry point; the scoring audit lives in src.scoring_validation.

Prefer ``make scoring-validation`` or ``python -m src.scoring_validation``.
The recorded ``python -m scripts.validate_scoring_activity`` command still works.
"""

from src.scoring_validation import main

if __name__ == "__main__":
    raise SystemExit(main())
