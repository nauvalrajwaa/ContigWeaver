"""Allow ``python -m contignexus`` invocation."""
from contignexus.pipeline import main
import sys

sys.exit(main())
