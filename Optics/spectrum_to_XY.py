#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  4 13:28:23 2026

@author: Igor Coropceanu
"""

import colour


def spectrum_to_XY(spectrum):
    """Compute CIE 1931 (x, y) chromaticity coordinates for a self-luminous spectrum.

    Parameters
    ----------
    spectrum : pandas.Series
        Spectral power distribution of an emissive source, indexed by
        wavelength (nm).

    Returns
    -------
    tuple[float, float]
        The (x, y) chromaticity coordinates.
    """
    sd = colour.SpectralDistribution(spectrum)

    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    illuminant = colour.SDS_ILLUMINANTS["E"]

    XYZ = colour.sd_to_XYZ(sd, cmfs, illuminant)
    x, y = colour.XYZ_to_xy(XYZ)

    return x, y
