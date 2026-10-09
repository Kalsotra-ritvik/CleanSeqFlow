"""Allow ``python -m transcriptomics_preprocessor``."""

import sys

from .cli import main

sys.exit(main())
