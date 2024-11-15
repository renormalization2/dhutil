# %%
from glob import glob
from tqdm import tqdm
import pandas as pd
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.io import fits
from astropy.table import Table
from astroquery.sdss import SDSS
from astropy import coordinates as coords
import astropy.units as u
import astropy.constants as const
import healpy as hp
import astropy_healpix as ah
from tqdm import tqdm
from matplotlib.patches import Circle
from astropy.coordinates import Angle
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy.cosmology.units as cu
from astropy.cosmology import (
    Planck18 as cosmo,
)  # Planck18 is a predefined ΛCDM cosmology
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon, Point


def select_skymap(num, mixed=False):
    datadir = Path("/Users/dhhyun/VSCode/GWInference/skymap")
    if not datadir.exists():
        # print("No", datadir)
        datadir = Path("/data3/dhhyun/GWInference/skymap")
    if not datadir.exists():
        # print("No", datadir)
        datadir = input("/data3/dhhyun/GWInference/skymap")
    O3 = "IGWN-GWTC3p0-v2-PESkyLocalizations"
    O2 = "IGWN-GWTC2p1-v2-PESkyMaps"

    O3_files = glob(str(datadir / O3 / f"*{num}*"))
    O2_files = glob(str(datadir / O2 / f"*{num}*"))
    globbed = [*O3_files, *O2_files]
    if len(globbed) > 1 and not mixed:
        strs = [Path(path).stem for path in globbed]
        popup = ""
        i = 0
        for s in strs:
            popup = popup + s + f"\t: {i}" + "\n"
            i += 1
        user_sel = int(input(popup))
        globbed = globbed[user_sel]
    elif mixed:
        globbed = globbed[1]
    else:
        globbed = globbed[0]
    return globbed


def plot_gw(event_num, ci=0.90, show=True, return_radec=False):
    # skymap = Table.read(datadir / file)
    skymap = Table.read(select_skymap(event_num))
    skymap.sort("PROBDENSITY", reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap["PROBDENSITY"]
    cumprob = np.cumsum(prob)
    # ci = 0.99  # 999  # 0.95
    i = cumprob.searchsorted(ci)
    # i = len(prob)
    skymap_ = skymap[:i]
    skymap_.sort("UNIQ")
    skymap_1 = skymap_["UNIQ",]
    uni = skymap_["UNIQ"]
    level, ipix = ah.uniq_to_level_ipix(uni)
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    # print(ra, dec)

    fig, ax = plt.subplots(dpi=1000)
    ax.scatter(
        ra * 180 / np.pi,
        dec * 180 / np.pi,
        label=f"up to {ci*100:.0f}%",
        marker=".",
        s=5,
    )
    stepsize = 5
    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(np.arange(int(start), end + stepsize, stepsize))
    stepsize = 5
    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(np.arange(int(start), end + stepsize, stepsize))
    ax.set_xlabel("RA ($^\circ$)")
    ax.set_ylabel("Dec ($^\circ$)")
    ax.set_aspect("equal")
    ax.legend()
    return fig, ax


def return_radec(event_num, ci=0.90, mixed=True):
    # skymap = Table.read(datadir / file)
    skymap = Table.read(select_skymap(event_num, mixed=mixed))  # Force 1
    skymap.sort("PROBDENSITY", reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap["PROBDENSITY"]
    cumprob = np.cumsum(prob)
    # ci = 0.99  # 999  # 0.95
    i = cumprob.searchsorted(ci)
    # i = len(prob)
    skymap_ = skymap[:i]
    skymap_.sort("UNIQ")
    skymap_1 = skymap_["UNIQ",]
    uni = skymap_["UNIQ"]
    level, ipix = ah.uniq_to_level_ipix(uni)
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    return ra, dec


def get_gw(event_num, ci=0.90, show=True, return_radec=False, use_mixed=True):
    """
    when return_radec is true, returns (ra , dec) in rad
    """
    # skymap = Table.read(datadir / file)
    skymap = Table.read(select_skymap(event_num, mixed=use_mixed))
    skymap.sort("PROBDENSITY", reverse=True)
    level, ipix = ah.uniq_to_level_ipix(skymap["UNIQ"])
    pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
    prob = pixel_area * skymap["PROBDENSITY"]
    cumprob = np.cumsum(prob)
    # ci = 0.99  # 999  # 0.95
    i = cumprob.searchsorted(ci)
    # i = len(prob)
    skymap_ = skymap[:i]
    skymap_.sort("UNIQ")
    skymap_1 = skymap_["UNIQ",]
    uni = skymap_["UNIQ"]
    level, ipix = ah.uniq_to_level_ipix(uni)
    nside = ah.level_to_nside(level)
    ra, dec = ah.healpix_to_lonlat(ipix, nside, order="nested")
    # print(ra, dec)

    if show:
        fig, ax = plt.subplots(dpi=1000)
        ax.scatter(
            ra * 180 / np.pi,
            dec * 180 / np.pi,
            label=f"up to {ci*100:.0f}%",
            marker=".",
            s=5,
        )
        stepsize = 5
        start, end = ax.get_xlim()
        ax.xaxis.set_ticks(np.arange(int(start), end + stepsize, stepsize))
        stepsize = 5
        start, end = ax.get_ylim()
        ax.yaxis.set_ticks(np.arange(int(start), end + stepsize, stepsize))
        ax.set_xlabel("RA ($^\circ$)")
        ax.set_ylabel("Dec ($^\circ$)")
        ax.set_aspect("equal")
        ax.legend()

    if return_radec:
        return ra, dec
    else:
        return fig, ax


class a_u:
    def __init__(self, nominal, perr, merr):
        merr = np.abs(merr)
        self.nominal = nominal
        self.upper = nominal + perr
        self.lower = nominal - merr
        self.array = np.array([nominal - merr, nominal, nominal + perr])


if __name__ == "__main__":
    plot_gw("190814")
