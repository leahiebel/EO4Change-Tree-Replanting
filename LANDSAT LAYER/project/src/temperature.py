"""Landsat Collection 2 surface-temperature conversion utilities."""

from __future__ import annotations

import ee

ST_B10_SCALE = 0.00341802
ST_B10_OFFSET_K = 149.0
KELVIN_TO_CELSIUS = 273.15


def st_b10_dn_to_celsius_value(digital_number: float) -> float:
    """Convert a Landsat ST_B10 digital number to degrees Celsius."""
    return digital_number * ST_B10_SCALE + ST_B10_OFFSET_K - KELVIN_TO_CELSIUS


def add_lst_celsius(image: ee.Image) -> ee.Image:
    """Add the scientific LST_C band from Landsat ST_B10.

    The output is not clamped to any display range. Visual min/max limits are
    applied only in web style configuration.
    """
    lst_c = (
        image.select("ST_B10")
        .multiply(ST_B10_SCALE)
        .add(ST_B10_OFFSET_K)
        .subtract(KELVIN_TO_CELSIUS)
        .rename("LST_C")
        .toFloat()
    )
    return image.addBands(lst_c, overwrite=True)

