"""Allow ``python -m contigweaver`` invocation."""
from contigweaver.pipeline import main
import sys

sys.exit(main())
