# dhutil
This is a collection of miscellaneous code snippets by Donghwan Hyun. 
The repository is public but is primarily intended for team members of GWUniverse at SNU.


Projected visualizations of a wide sky area can be made more convenient with `ligo.skymap`, which can be installed as follows.
Also refer to their official documentation: https://lscsoft.docs.ligo.org/ligo.skymap/quickstart/install.html
```
$ conda config --add channels conda-forge
$ conda config --set channel_priority strict
$ conda install ligo.skymap
```



##  Install

Clone the repository

```
$ git clone https://github.com/renormalization2/dhutil.git
```

navigate into the cloned directory

```
$ cd dhutil/
```

and install it in editable mode.

```
$ pip install -e .
```

Now you can use it just like other packages. e.g.,
```
import matplotlib.pyplot as plt
from astropy.wcs import WCS
import dhutil as dh

# Define WCS Projection
wcs = WCS(naxis=2)
wcs.wcs.crval = [115, -27.0]  # Center of the projection in RA/Dec (degrees)
wcs.wcs.cdelt = [-0.01, 0.01]  # Pixel scale in degrees/pixel
wcs.wcs.crpix = [100, 100]  # Reference pixel (center of the plot)
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]  # Gnomonic (tangent plane) projection

# Visualize 7DS Tiles with a one-liner!
fig, ax = plt.subplots(subplot_kw={"projection": wcs}, dpi=300)
ax.plot([110, 120], [-25, -30], "o--", transform=ax.get_transform("world"))
ax.coords[0].set_format_unit('deg')  # convert ra unit (hour to deg)
ax.grid()

dh.set_xylim(ax, 108, 123, -34, -21)  # optional
dh.overlay_tiles(color="k", fontsize=6.6)
```


## Setting Up Paths to Large Data

If you use this code in QSO, the SNU server for astronomy, the absolute paths to the GW skymaps are already set.
But you can always look into `gwloc.select_skymap` if you are in another environment or want to add other skymaps.