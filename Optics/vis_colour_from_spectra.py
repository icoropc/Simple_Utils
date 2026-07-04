#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 10:00:10 2026

@author: Igor Coropceanu


"""

import numpy as np
import matplotlib.pyplot as plt
import colour

# ---------------------------------------------------------------------------
# Gaussian component definitions: (relative amplitude, FWHM [nm], peak [nm])
# ---------------------------------------------------------------------------
components = [
    {"amplitude": 1.0, "fwhm": 4.0, "peak": 462.0},  # blue
    {"amplitude": 0.000, "fwhm": 30.0, "peak": 540.0},  # green
    {"amplitude": 0.000, "fwhm": 30.0, "peak": 630.0},  # red
]

SHAPE = colour.SpectralShape(380, 780, 1)
wavelengths = SHAPE.wavelengths


def gaussian(wavelengths, amplitude, fwhm, peak):
    """Gaussian defined by peak amplitude and full width at half maximum."""
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return amplitude * np.exp(-0.5 * ((wavelengths - peak) / sigma) ** 2)


def build_spectrum(components, wavelengths):
    """Return the individual component arrays and their summed spectrum."""
    component_values = [
        gaussian(wavelengths, c["amplitude"], c["fwhm"], c["peak"])
        for c in components
    ]
    total = np.sum(component_values, axis=0)
    return component_values, total


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


def plot_result(wavelengths, component_values, total, rgb, rgb_clipped):
    fig, (ax_spectrum, ax_swatch) = plt.subplots(
        1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [3, 1]}
    )

    colors = ["tab:blue", "tab:green", "tab:red"]
    for values, comp, color in zip(component_values, components, colors):
        ax_spectrum.plot(
            wavelengths,
            values,
            "--",
            color=color,
            linewidth=1.5
        )

    ax_spectrum.plot(
        wavelengths, total, color="black", linewidth=2.5, label="Total spectrum"
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
    component_values, total = build_spectrum(components, wavelengths)
    rgb, rgb_clipped, sd = spectrum_to_rgb(total, wavelengths)


    print(f"\nsRGB (raw):     {rgb}")
    print(f"sRGB (clipped): {rgb_clipped}")
    print(f"Hex:            {colour.notation.RGB_to_HEX(rgb_clipped)}")

    plot_result(wavelengths, component_values, total, rgb, rgb_clipped)


if __name__ == "__main__":
    main()
