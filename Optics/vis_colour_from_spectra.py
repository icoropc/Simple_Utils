#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 10:00:10 2026

@author: Igor Coropceanu


"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import colour

# ---------------------------------------------------------------------------
# Gaussian component definitions: (relative amplitude, FWHM [nm], peak [nm])
# ---------------------------------------------------------------------------
components = {
    "blue_spectrum": {"amplitude": 1.0, "fwhm": 4.0, "peak": 457.0},
    "green_spectrum": {"amplitude": 0.6, "fwhm": 40.0, "peak": 540.0},
    "red_spectrum": {"amplitude": 0.6, "fwhm": 40.0, "peak": 630.0},
}

SHAPE = colour.SpectralShape(380, 780, 1)
wavelengths = SHAPE.wavelengths


def gaussian(wavelengths, amplitude, fwhm, peak):
    """Gaussian defined by peak amplitude and full width at half maximum."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((wavelengths - peak) / sigma) ** 2)


def build_spectrum(components, wavelengths):
    """Return each component as a named pandas Series and their summed spectrum."""
    spectra = {
        name: pd.Series(
            gaussian(wavelengths, c["amplitude"], c["fwhm"], c["peak"]),
            index=wavelengths,
            name=name,
        )
        for name, c in components.items()
    }
    total_spectrum = sum(spectra.values())
    total_spectrum.name = "total_spectrum"
    return spectra, total_spectrum


def spectrum_to_rgb(values, wavelengths, illuminant_sd=None):
    """Convert a self-luminous spectral power distribution to sRGB (0-1 clipped).

    The modelled spectrum is emitted (e.g. an LED), not reflected, so it
    must not be re-weighted by a real-world illuminant like D65 - that
    would treat the LED as a reflective sample lit by a second light
    source. The default illuminant is therefore CIE Illuminant E
    (equal-energy), which leaves the spectrum's shape untouched.
    """
    sd = colour.SpectralDistribution(
        dict(zip(wavelengths, values)), name="Modelled spectrum"
    )

    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"]
    illuminant = (
        illuminant_sd
        if illuminant_sd is not None
        else colour.SDS_ILLUMINANTS["E"]
    )

    XYZ = colour.sd_to_XYZ(sd, cmfs, illuminant) / 100.0
    illuminant_xy = colour.CCS_ILLUMINANTS[
        "CIE 1931 2 Degree Standard Observer"
    ]["E"]

    RGB = colour.XYZ_to_sRGB(XYZ, illuminant_xy)
    RGB_normalized = np.clip(RGB / np.max(RGB), 0.0, 1.0)

    return RGB, RGB_normalized, sd


def plot_result(wavelengths, spectra, total_spectrum, rgb, rgb_clipped):
    fig, (ax_spectrum, ax_swatch) = plt.subplots(
        1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [3, 1]}
    )

    for name, spectrum in spectra.items():
        color = f"tab:{name.split('_')[0]}"
        ax_spectrum.plot(
            spectrum,
            "--",
            color=color,
            linewidth=1.5
        )

    ax_spectrum.plot(
        total_spectrum, color="black", linewidth=2.5, label="Total spectrum"
    )
    ax_spectrum.set_xlabel("Wavelength (nm)")
    ax_spectrum.set_ylabel("Relative power")
    ax_spectrum.set_title("Modelled spectrum (sum of 3 Gaussians)")
    ax_spectrum.legend(fontsize=8, loc="upper right")
    ax_spectrum.set_xlim(wavelengths.min(), wavelengths.max())
    ax_spectrum.grid(alpha=0.3)

    ax_swatch.imshow([[rgb_clipped]])
    ax_swatch.set_xticks([])
    ax_swatch.set_yticks([])
    hex_code = colour.notation.RGB_to_HEX(rgb_clipped)
    ax_swatch.set_title(
        "Resulting sRGB colour\n"
        f"RGB (raw): ({rgb[0]:.3f}, {rgb[1]:.3f}, {rgb[2]:.3f})\n"
        f"RGB (clipped): ({rgb_clipped[0]:.3f}, {rgb_clipped[1]:.3f}, "
        f"{rgb_clipped[2]:.3f})\n"
        f"Hex: {hex_code}",
        fontsize=9,
    )

    fig.tight_layout()
    plt.show()


def main():
    spectra, total_spectrum = build_spectrum(components, wavelengths)
    red_spectrum = spectra["red_spectrum"]
    green_spectrum = spectra["green_spectrum"]
    blue_spectrum = spectra["blue_spectrum"]

    rgb, rgb_clipped, sd = spectrum_to_rgb(total_spectrum, wavelengths)

    print(f"\nsRGB (raw):     {rgb}")
    print(f"sRGB (clipped): {rgb_clipped}")
    print(f"Hex:            {colour.notation.RGB_to_HEX(rgb_clipped)}")

    plot_result(wavelengths, spectra, total_spectrum, rgb, rgb_clipped)


if __name__ == "__main__":
    main()
