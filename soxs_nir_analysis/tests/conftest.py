import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

# MAKE THE SCRIPT IMPORTABLE AS A MODULE.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
