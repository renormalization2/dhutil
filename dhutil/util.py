# %%


def h2deg(h):
    import astropy.units as u

    deg = (h * u.hourangle).to(u.deg).value
    return deg


def deg2h(deg):
    import astropy.units as u

    h = (deg * u.deg).to(u.hourangle).value
    return h


def deg2hms(deg):
    from astropy.coordinates import Angle

    ang = Angle(deg, unit="deg")
    hms = ang.hms
    # h = (deg * u.deg).to(u.hourangle).value
    return hms.h, hms.m, round(hms.s, 6)


# def lapse(cp1=None, explanation="elapsed", report=True, end="\n"):
#     """
#     cp for CheckPoint!
#     lapse takes a predefined timeit.default_timer value
#     and returns the current timeit.default_timer value,
#     printing the difference of the two with an explanation.

#     use it like
#     cp = dh.lapse()
#     cp = dh.lapse(cp, 'elapsed for your task')
#     """
#     from timeit import default_timer as timer

#     if cp1 == None:
#         cp1 = timer()
#         return cp1

#     cp2 = timer()
#     if report:
#         dt = cp2 - cp1
#         # dt = dt if dt < 60 else (dt / 3600 if dt > 3600 else dt / 60)
#         if dt < 60:
#             dt, u = dt, "s"
#         elif dt > 3600:
#             dt, u = dt / 3600, "h"
#         else:
#             dt, u = dt / 60, "m"
#         print(f"{dt:.3f}{u} {explanation}", end=end)
#     return cp2


# Global variable for checkpoint storage
_dhutil_lapse_checkpoint = None


def lapse(explanation="elapsed", report=True, end="\n"):
    """
    A utility function to measure and report elapsed time using a global checkpoint.

    Parameters:
    explanation (str): Description for the elapsed time report.
    report (bool): Whether to print the elapsed time.
    end (str): Ending character for the printed output.

    Usage:
    lapse("Start")  # Initializes the timer and prints "Timer started"
    lapse("Task completed")  # Prints the time elapsed since the last call
    """
    from timeit import default_timer as timer

    global _dhutil_lapse_checkpoint  # Access or modify the global variable

    current_time = timer()

    if _dhutil_lapse_checkpoint is None:  # Initialize if it's the first call
        _dhutil_lapse_checkpoint = current_time
        # if report:
        #     print(f"Timer started: {explanation}", end=end)
    else:
        elapsed_time = current_time - _dhutil_lapse_checkpoint

        if elapsed_time < 60:
            dt, unit = elapsed_time, "s"
        elif elapsed_time > 3600:
            dt, unit = elapsed_time / 3600, "h"
        else:
            dt, unit = elapsed_time / 60, "m"
        if report:
            print(f"{dt:.3f}{unit} {explanation}", end=end)
        _dhutil_lapse_checkpoint = current_time  # Update the checkpoint
        return dt


# class Lapse:
#     """
#     class ver. of the function lapse
#     """

#     def __init__(self):
#         from timeit import default_timer as timer

#         self.init_time = timer()

#     def update(self):
#         if cp1 == None:
#             cp1 = timer()
#             return cp1

#         cp2 = timer()
#         if report:
#             dt = cp2 - cp1
#             dt = dt if dt < 60 else (dt / 3600 if dt > 3600 else dt / 60)
#             print(f" {dt:.3f}s {explanation}", end=end)
#         return cp2


def loghist(arr, binsize=30, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt

    _min, _max = np.min(arr), np.max(arr)
    if _min >= 0:
        bins = np.logspace(np.log10(_min), np.log10(_max), binsize)
        plt.hist(arr, bins=bins, **kwargs)
    else:
        raise ValueError


def rangecut(column, lower_cut, upper_cut):
    condition = (lower_cut < column) & (column < upper_cut)
    return condition
    # return table[condition]


def save(var, filename="data.pkl"):
    import pickle

    with open(filename, "wb") as f:
        pickle.dump(var, f)


def load(filename="data.pkl"):
    import pickle

    with open(filename, "rb") as f:
        return pickle.load(f)


def read(filename="data.pkl"):
    load(filename)


def mjd2moonphase(mjd):
    import numpy as np
    from astropy.time import Time
    from astropy.coordinates import get_sun, get_body

    if isinstance(mjd, Time):
        time = mjd
    elif isinstance(mjd, str):
        # Convert MJD to an astropy Time object
        time = Time(mjd, format="mjd")

    # Get positions of the Sun and Moon at this time
    sun_pos = get_sun(time)
    moon_pos = get_body("moon", time)

    # Calculate the angle between the Sun and Moon in degrees
    phase_angle = sun_pos.separation(moon_pos).deg

    # Calculate Moon phase as a percentage based on the phase angle
    phase = (1 + np.cos(np.radians(phase_angle))) / 2 * 100  # Illuminated fraction in %

    # Output the moon phase as a percentage
    return phase


def viewer_data(img, dpi=100, cmap="viridis", save=None, **kwargs):
    import numpy as np
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from astropy.visualization import ZScaleInterval

    fig, ax = plt.subplots(1, dpi=dpi)  # , figsize=(5,5))
    if "lim" in kwargs:
        xl, xu, yl, yu = kwargs["lim"]
    else:
        xl, xu, yl, yu = 0, np.shape(img)[1], 0, np.shape(img)[0]
    img = img[yl:yu, xl:xu]

    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(img)
    im = ax.imshow(img, vmin=vmin, vmax=vmax, cmap=cmap)  # origin='lower',
    # ax.set_title(title[i])
    ax.tick_params(axis="both", length=0.0, labelleft=False, labelbottom=False)
    # ax.text(0.05, 0.95, label[i], fontsize=15.0, fontweight='bold',
    #         transform=ax.transAxes, ha='left', va='top')

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    # plt.colorbar(im, ax=ax)
    # plt.show()
    if save is not None:
        plt.imsave(save, np.clip(img, 1.03 * vmin, 0.97 * vmax), cmap=cmap)


def unique(base_filename):
    """
    Modify base_filename by appending a number to make it unique.
    If base_filename is 'image.png' and it exists, this function will return
    'image_1.png', 'image_2.png', etc.
    """
    import os

    counter = 1
    filename, file_extension = os.path.splitext(base_filename)
    new_filename = base_filename

    while os.path.exists(new_filename):
        new_filename = f"{filename}_{counter}{file_extension}"
        print(new_filename, "already exists")
        counter += 1

    print("saving as", new_filename)
    return new_filename


def str2bool(val):
    """Convert a string representation of truth to true (1) or false (0).
    True values are 'y', 'yes', 't', 'true', 'on', and '1'; false values
    are 'n', 'no', 'f', 'false', 'off', and '0'.  Raises ValueError if
    'val' is anything else.
    """
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return 1
    elif val in ("n", "no", "f", "false", "off", "0"):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))


# def save_as_fits_ldac(outbl, tablename):
#     """
#     Save input table as FITS_LDAC format

#     Parameters:
#     outbl (Table): Astropy Table to save
#     tablename (str): save name

#     Returns:
#     None
#     """
#     from astropy.io import fits

#     # BINTABLE
#     hdul = fits.HDUList()

#     # Initialize Primary HDU (Header Only, no data)
#     primary_hdu = fits.PrimaryHDU()
#     hdul.append(primary_hdu)

#     # input table to LDAC_OBJECTS
#     bintable_hdu = fits.BinTableHDU(outbl, name="LDAC_OBJECTS")
#     hdul.append(bintable_hdu)

#     # save as fits ldac
#     hdul.writeto(tablename, overwrite=True)


def write_ldac(catalog_table, filename):
    """
    Save an astropy Table as a valid FITS-LDAC catalog
    that SCAMP / SWarp / etc. will accept.
    """
    from astropy.io import fits

    # 1) Primary HDU (empty)
    hdu_primary = fits.PrimaryHDU()

    # 2) LDAC_IMHEAD extension (contains image header cards)
    #    SExtractor typically writes one row per 80-char "card".
    #    Here we make a trivial/empty version.
    col_imhead = fits.Column(
        name="Field Header Card", format="80A", array=[b""] * 1
    )  # At least 1 row
    hdu_imhead = fits.BinTableHDU.from_columns([col_imhead], name="LDAC_IMHEAD")

    # 3) LDAC_OBJECTS extension (your actual catalog data)
    hdu_objects = fits.BinTableHDU(catalog_table, name="LDAC_OBJECTS")

    # Combine them into an HDUList
    hdulist = fits.HDUList([hdu_primary, hdu_imhead, hdu_objects])

    # Write to disk
    hdulist.writeto(filename, overwrite=True)


def overlay_hptiles(nside=32, ra_min=30, ra_max=40, dec_min=-40, dec_max=-30):
    import numpy as np
    import matplotlib.pyplot as plt
    import healpy as hp

    # Total number of pixels in this resolution
    npix = hp.nside2npix(nside)

    # Plot setup
    ax = plt.gca()
    ax.set_xlabel("Right Ascension (deg)")
    ax.set_ylabel("Declination (deg)")
    # ax.grid(True)

    # Loop over all HEALPix pixels to plot their edges
    for pix in range(npix):
        # Get the boundary vertices of the pixel in 3D Cartesian coordinates
        boundaries = hp.boundaries(nside, pix, step=1)  # shape (3, n_vertices)

        # Convert 3D Cartesian coordinates to spherical coordinates (theta, phi)
        x, y, z = boundaries
        theta = np.arccos(z)  # theta = arccos(z)
        phi = np.arctan2(y, x)  # phi = atan2(y, x)

        # Convert to RA, DEC (in degrees)
        ra = np.rad2deg(phi)  # phi corresponds to RA
        dec = np.rad2deg(0.5 * np.pi - theta)  # theta to DEC

        # Ensure RA is within [0, 360] degrees
        # error handling for ra ~ 0 to be added later
        ra = (ra + 360) % 360

        # Check if the tile is within the desired RA and DEC range
        in_ra = ra.min() >= ra_min and ra.max() <= ra_max
        in_dec = dec.min() >= dec_min and dec.max() <= dec_max
        if in_ra and in_dec:
            # Plot the edges as a polygon (close the loop by appending the first point at the end)
            plt.plot(
                np.append(ra, ra[0]),
                np.append(dec, dec[0]),
                color="blue",
                linewidth=0.5,
            )
            # Annotate the tile with the pixel index
            plt.text(
                ra.mean(),
                dec.mean(),
                str(pix),
                fontsize=6,
                color="red",
                ha="center",
                va="center",
            )


# %%
if __name__ == "__main__":
    lapse()
    print("a")
    lapse()
