from typing import Iterable, Tuple, Sequence, Any
import numpy as np
from astropy.coordinates import SkyCoord, Distance
from astropy.time import Time
from astropy.table import Table, hstack, MaskedColumn
import astropy.units as u
from astropy.units import Quantity
import operator


def match_two_catalogs(
    sci_tbl: Table,
    ref_tbl: Table,
    *,
    x0: str = "ALPHA_J2000",
    y0: str = "DELTA_J2000",
    x1: str | None = None,
    y1: str | None = None,
    radius: float | Quantity = 1,
    how: str = "inner",
    correct_pm: bool = False,
    obs_time: Time | None = None,
    pm_keys: dict | None = None,
) -> Table:
    """
    Cross-match two catalogues on the sky and (optionally) apply proper-motion
    correction to the *reference* catalogue before matching.
    The coordinates are assumed to be Equatorial (RA, Dec) in degrees.

    Parameters
    ----------
    sci_tbl, ref_tbl
        `astropy.table.Table` objects to be matched.
    x0, y0, x1, y1
        Column names that contain right ascension and declination *in degrees*.
        If `x1` / `y1` are omitted they default to `x0` / `y0`.
    radius
        Maximum separation for a match.  A bare ``float`` is interpreted
        in **arcseconds**; a `~astropy.units.Quantity` may carry any angle unit.
    how : {'inner', 'left'}, optional
        Join strategy:

        * ``'inner'`` - return only matched rows (default)
        * ``'left'``  - return every row of *sci_tbl* and mask unmatched
                        entries from ref_tbl
    correct_pm
        If *True*, propagate stars in the *reference* catalogue from their
        catalogued epoch (`pm_info['ref_epoch']`, default 2016.0 TDB) to
        `obs_time` using their proper motion and parallax.
    obs_time
        Observation time of the *science* catalogue as an `~astropy.time.Time`
        instance.  Required when ``correct_pm=True``.
    pm_info
        Mapping that defines column names for the kinematic quantities and the
        reference epoch::

            dict(pmra="pmra", pmdec="pmdec",
                 parallax="parallax", ref_epoch=2016.0)

        You can pass additional keys (e.g. ``"rv"``) without hurting anything.

    Returns
    -------
    Table
        A merged `~astropy.table.Table` with one row per successful positional
        match.  Columns from the smaller (“driver”) catalogue keep their
        original names; duplicate names coming from the second catalogue receive
        the suffix ``"_ref"`` or ``"_sci"`` depending on which side they came
        from.

    Notes
    -----
    * Rows whose kinematic information is incomplete are **kept** but **not**
      epoch-propagated.  They will still match if their uncorrected position
      lies within *radius*.
    * Negative or zero parallaxes are converted into ``NaN`` distances via
      ``Distance(parallax, allow_negative=True)``; these rows are likewise
      retained but un-moved.
    """

    if how not in {"inner", "left"}:
        raise ValueError("how must be 'inner' or 'left'")

    if x1 is None:
        x1 = x0
    if y1 is None:
        y1 = y0

    coord_sci = SkyCoord(sci_tbl[x0], sci_tbl[y0], unit="deg", copy=False)
    coord_ref = SkyCoord(ref_tbl[x1], ref_tbl[y1], unit="deg", copy=False)

    if correct_pm:
        from astropy.coordinates import Distance

        if obs_time is None:
            raise ValueError("obs_time must be provided if correct_pm is True")

        if pm_keys is None:
            pm_keys = dict(pmra="pmra", pmdec="pmdec", parallax="parallax", ref_epoch=2016.0)

        # Vectorised columns with units
        pm_ra = ref_tbl[pm_keys["pmra"]] * u.mas / u.yr
        pm_dec = ref_tbl[pm_keys["pmdec"]] * u.mas / u.yr
        dist = Distance(parallax=ref_tbl[pm_keys["parallax"]] * u.mas, allow_negative=True)  # allow (-) parallax

        good = np.isfinite(pm_ra) & np.isfinite(pm_dec) & np.isfinite(dist)

        if np.any(good):  # only build a 6-D frame when possible
            moved = SkyCoord(
                ra=ref_tbl[x1][good] * u.deg,
                dec=ref_tbl[y1][good] * u.deg,
                pm_ra_cosdec=pm_ra[good],
                pm_dec=pm_dec[good],
                distance=dist[good],
                obstime=Time(pm_keys["ref_epoch"], format="jyear"),
            ).apply_space_motion(new_obstime=obs_time)

            # update with the moved coordinates
            # coord_ref = coord_ref.replicate(copy=False, obstime=obs_time)  # this reuses memory buffer and is efficient
            coord_ref.ra[good] = moved.ra
            coord_ref.dec[good] = moved.dec
            coord_ref._sky_coord_frame.cache.clear()

    # ---------------- driver selection --------------------------------
    if how == "left":
        # science catalogue must drive so we can return all its rows
        coord0, coord1 = coord_sci, coord_ref
        tbl0, tbl1 = sci_tbl, ref_tbl
        tag1 = "ref"
    elif len(coord_sci) > len(coord_ref):  # small-driver heuristic for a simple common table
        coord0, coord1 = coord_ref, coord_sci  # match_to_catalog_sky is N0 * logN1 complexity
        tbl0, tbl1 = ref_tbl, sci_tbl
        tag1 = "sci"
    else:
        coord0, coord1 = coord_sci, coord_ref
        tbl0, tbl1 = sci_tbl, ref_tbl
        tag1 = "ref"

    # --- catalogue matching & merge ------------------------------------------
    idx, sep2d, _ = coord0.match_to_catalog_sky(coord1)
    rtol = radius * u.arcsec if not isinstance(radius, Quantity) else radius
    matched = sep2d < rtol

    # This works too...
    # merged = sci_tbl.join(ref_tbl, keys_left=coord_sci, keys_right=coord_ref,
    #                   rtol=0*u.arcsec + radius, table_names=["sci", "ref"])

    # ------------------------------- keep all science rows -------------
    if how == "left":
        out = tbl0.copy()  # shallow copy; cheap
        ref_slice = tbl1[idx]  # reference rows in match order

        # disambiguate duplicate column names
        dupes = set(out.colnames) & set(ref_slice.colnames)
        for name in dupes:
            ref_slice.rename_column(name, f"{name}_{tag1}")

        # add reference columns as masked arrays where no match
        for name in ref_slice.colnames:
            out[name] = MaskedColumn(ref_slice[name].data, mask=~matched)

        # separation column (arcsec)
        out["separation"] = MaskedColumn(sep2d.arcsec, mask=~matched)
        return out

    # ------------------------------- default: return only matches ------
    else:
        m0 = tbl0[matched]
        m1 = tbl1[idx[matched]]

        dupes = set(m0.colnames) & set(m1.colnames)
        for name in dupes:
            m1.rename_column(name, f"{name}_{tag1}")

        merged = hstack([m0, m1], join_type="exact")
        merged["separation"] = sep2d[matched].arcsec
        return merged


Condition = Tuple[str, Any, str]  # (column, value, method)

_OPS = {
    # greater-than
    "lower": operator.gt,  ">":  operator.gt,
    # greater-or-equal
    ">=":   operator.ge,
    # less-than
    "upper": operator.lt,  "<":  operator.lt,
    # less-or-equal
    "<=":   operator.le,
    # equal
    "equal": operator.eq,  "==": operator.eq,  "=": operator.eq,
}  # fmt: skip


def build_condition_mask(table, conditions: Iterable[Condition]) -> np.ndarray:
    """
    Return a boolean mask that is True only for rows satisfying *all* conditions.

    Parameters
    ----------
    table : Table | DataFrame | structured ndarray
        Object supporting ``table[col]`` column access.
    conditions : iterable of (key, method, value)
        method can be any alias listed in _METHOD_MAP (case-insensitive).

    Returns
    -------
    numpy.ndarray[bool]  shape (len(table),)
    """

    # later incorporate numexpr

    conditions = _parse_conditions(conditions)
    mask = np.ones(len(table), dtype=bool)

    for key, method, value in conditions:
        m = method.strip().lower()
        if m not in _OPS:
            raise ValueError(f"Unknown method '{method}'. Allowed: {', '.join(_OPS)}")
        mask &= _OPS[m](table[key], value)

    return mask


def _parse_conditions(conditions: Sequence[Any]) -> Iterable[Condition]:
    """
    Turn *raw* into an iterable of (key, op, value) tuples.

    *raw* may be:
      • already an iterable of 3-tuples
      • a flat 1-D list/array whose length is a multiple of 3
    """
    # Understand something like [(k, op, v), …]
    if conditions and isinstance(conditions[0], (tuple, list)) and len(conditions[0]) == 3:
        return conditions  # type: ignore[arg-type]

    # If not, try the “flat” form
    if len(conditions) % 3 != 0:
        raise ValueError("Flat conditions must have length divisible by 3 " "(key, op, value repeated).")
    it = iter(conditions)
    return list(zip(it, it, it))


def filter_table(table: Table, conditions: Iterable[Condition]) -> Table:
    mask = build_condition_mask(table, conditions)
    return table[mask]
