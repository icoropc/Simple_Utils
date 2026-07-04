#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 13:28:23 2026

@author: Igor Coropceanu
"""

import numpy as np
import pandas as pd
import colour

DEFAULT_SHAPE = colour.SpectralShape(380, 780, 1)


def spectrum_to_XY(spectrum):
    """Compute CIE 1931 (x, y) chromaticity coordinates for a self-luminous spectrum.

    Parameters
    ----------
    spectrum : pandas.Series or numpy.ndarray
        Spectral power distribution of an emissive source. If a pandas
        Series, its index is used as the wavelengths (nm); if a numpy
        array, wavelengths default to DEFAULT_SHAPE (380-780nm, 1nm step).

    Returns
    -------
    tuple[float, float]
        The (x, y) chromaticity coordinates.
    """
    if isinstance(spectrum, pd.Series):
        wavelengths = spectrum.index.to_numpy()
        values = spectrum.to_numpy()
    else:
        values = np.asarray(spectrum)
        wavelengths = DEFAULT_SHAPE.wavelengths

    sd = colour.SpectralDistribution(dict(zip(wavelengths, values)))

    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    illuminant = colour.SDS_ILLUMINANTS["E"]

    XYZ = colour.sd_to_XYZ(sd, cmfs, illuminant)
    x, y = colour.XYZ_to_xy(XYZ)

    return x, y
