from .util import *
from .SDS import get_tiles, overlay_tiles, set_xylim

try:
    from .gwloc import *

    # from .gwloc import select_skymap, plot_gw, return_radec
    # from .gwloc import a_u
    from .gwloc import *  # get_gw
except ImportError as e:
    print(e)
    print(
        "Warning: 'gwloc' features are unavailable because the dependency is not installed."
    )
    print("This shouldn't matter as long as you don't use gwloc")

# from .SDS import plot_grid  # deprecated


# class PathManager:
#     def __init__(self):
#         from pathlib import Path
#         paths_to_try = [["abc", "abc"], ["abc", "abc"]]
