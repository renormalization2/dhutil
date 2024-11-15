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

    ang = Angle(deg, unit=u.deg)
    hms = ang.hms
    # h = (deg * u.deg).to(u.hourangle).value
    return hms.h, hms.m, round(hms.s, 6)


def lapse(cp1=None, explanation="elapsed", report=True, end="\n"):
    """
    lapse takes a predefined timeit.default_timer value
    and returns the current timeit.default_timer value,
    printing the difference of the two with an explanation.

    use it like
    cp = dh.lapse()
    cp = dh.lapse(cp, 'elapsed for your task')
    """
    from timeit import default_timer as timer

    if cp1 == None:
        cp1 = timer()
        return cp1

    cp2 = timer()
    if report:
        dt = cp2 - cp1
        dt = dt if dt < 60 else (dt / 3600 if dt > 3600 else dt / 60)
        print(f" {dt:.3f}s {explanation}", end=end)
    return cp2


class Lapse:
    """
    class ver. of the function lapse
    """

    def __init__(self):
        from timeit import default_timer as timer

        self.init_time = timer()

    def update(self):
        if cp1 == None:
            cp1 = timer()
            return cp1

        cp2 = timer()
        if report:
            dt = cp2 - cp1
            dt = dt if dt < 60 else (dt / 3600 if dt > 3600 else dt / 60)
            print(f" {dt:.3f}s {explanation}", end=end)
        return cp2


def save():
    import pickle


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
