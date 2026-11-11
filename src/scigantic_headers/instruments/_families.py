"""Imports every instrument family module so its sniffers register. Keep the
list sorted; each family adds one line."""

from . import cells  # noqa: F401
from . import flow_map_lh  # noqa: F401
from . import pcr  # noqa: F401
from . import plate_b  # noqa: F401
from . import spectro  # noqa: F401
