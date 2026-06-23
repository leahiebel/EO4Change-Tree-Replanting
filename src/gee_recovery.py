"""
EO4Change Group 4 — forest mask, RRI, recovery_class, guard_flags, sensitivity.
"""
import ee
from gee_config import (aoi, exp, FM_ENABLED, FM_YEAR, FM_MONTHS, FM_THRESHOLD,
                         T_RRI_GOOD, T_RRI_LOW, T_RRI_MIN_DENOM, T_BURN,
                         T_NDRE_GOOD, T_DNDVI_GOOD, T_BSI_HIGH, T_WATER, T_DISTURB)

# ── Class metadata ────────────────────────────────────────────────────────────
RECOVERY_CLASSES = [
    (0, "Uncertain / water / insufficient data", "#888888"),
    (1, "Recovering well",                       "#1f5e1f"),
    (2, "Recovering, but weak",                 "#a8d666"),
    (3, "Not recovering / failed",              "#ff1a1a"),
    (4, "Outside burn recovery assessment",     "#4b2e1f"),
]
RECOVERY_PALETTE = [c for _, _, c in sorted(RECOVERY_CLASSES)]
VIS_RECOVERY_CLASS = {"min": 0, "max": 4, "palette": RECOVERY_PALETTE}

GUARD_CLASSES = [
    (0, "Clear recovery signal",              "#00b894"),  # teal-green  — good signal
    (1, "Water / insufficient observations",  "#636e72"),  # neutral grey — no data
    (2, "Possible later disturbance",         "#d63031"),  # vivid red    — warning
    (3, "Mixed recovery signal",              "#fdcb6e"),  # amber        — caution
    (4, "Outside burn recovery assessment",   "#b2bec3"),  # light grey   — excluded
]
GUARD_PALETTE = [c for _, _, c in sorted(GUARD_CLASSES)]
VIS_GUARDS = {"min": 0, "max": 4, "palette": GUARD_PALETTE}


def build_forest_mask() -> "ee.Image | None":
    if not FM_ENABLED:
        return None
    ws = f"{FM_YEAR}-{FM_MONTHS[0]:02d}-01"
    we = f"{FM_YEAR}-{FM_MONTHS[1]:02d}-01"
    dw = (ee.ImageCollection("GOOGLE/DYNAMICWORLD/V1")
          .filterBounds(aoi)
          .filterDate(ws, we)
          .select("trees")
          .median()
          .clip(aoi))
    return dw.gte(FM_THRESHOLD).rename("forest_mask")


def _recovery_valid_mask(recent: ee.Image, change: ee.Image) -> ee.Image:
    return (recent.select("NDVI").mask()
            .multiply(recent.select("NDRE").mask())
            .multiply(recent.select("BSI").mask())
            .multiply(recent.select("NDWI").mask())
            .multiply(change.select("NDVI_change").mask())
            .multiply(change.select("NBR_change").mask()))


def compute_rri(recent: ee.Image, postfire: ee.Image,
                prefire: ee.Image) -> ee.Image:
    ndvi_recent   = recent.select("NDVI")
    ndvi_postfire = postfire.select("NDVI")
    ndvi_prefire  = prefire.select("NDVI")
    denom = ndvi_prefire.subtract(ndvi_postfire)
    rri = (ndvi_recent.subtract(ndvi_postfire)
                      .divide(denom.max(ee.Image.constant(T_RRI_MIN_DENOM)))
                      .clamp(-1, 2)
                      .rename("RRI")
                      .updateMask(denom.gte(T_RRI_MIN_DENOM)))
    return rri.clip(aoi)


def build_recovery_class(recent: ee.Image, postfire: ee.Image,
                         prefire: ee.Image, change: ee.Image) -> ee.Image:
    rri = compute_rri(recent, postfire, prefire)
    dnbr  = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(T_BURN)
    ndre_good = recent.select("NDRE").gte(T_NDRE_GOOD)
    water     = recent.select("NDWI").gt(T_WATER)
    valid     = _recovery_valid_mask(recent, change)
    good   = rri.gte(T_RRI_GOOD).And(ndre_good).And(water.Not()).And(burned)
    failed = rri.lt(T_RRI_LOW).And(water.Not()).And(burned)
    cls = ee.Image.constant(2).rename("recovery_class")
    cls = cls.where(good, 1)
    cls = cls.where(failed, 3)
    cls = cls.where(valid.eq(1).And(water.Not()).And(burned.Not()), 4)
    cls = cls.where(valid.eq(0).Or(water), 0)
    return cls.clip(aoi).toUint8().rename("recovery_class")


def build_guard_flags(recent: ee.Image, postfire: ee.Image,
                      prefire: ee.Image, change: ee.Image) -> ee.Image:
    dnbr = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(T_BURN)
    water = recent.select("NDWI").gt(T_WATER)
    valid = _recovery_valid_mask(recent, change)
    disturbance = (
        change.select("NBR_change")
        .lt(-T_DISTURB)
        .And(burned)
        .And(water.Not())
        .And(valid.eq(1))
    )
    mixed_signal = (
        change.select("NDVI_change").gte(T_DNDVI_GOOD)
        .And(
            recent.select("NDRE").lt(T_NDRE_GOOD)
            .Or(recent.select("BSI").gt(T_BSI_HIGH))
        )
        .And(burned)
        .And(water.Not())
        .And(valid.eq(1))
    )
    out = ee.Image.constant(0).rename("guard_flags")
    out = out.where(valid.eq(1).And(water.Not()).And(burned.Not()), 4)
    out = out.where(disturbance, 2)
    out = out.where(mixed_signal.And(disturbance.Not()), 3)
    out = out.where(valid.eq(0).Or(water), 1)
    return out.clip(aoi).toUint8().rename("guard_flags")


def build_recovery_class_sensitivity(
    recent: ee.Image,
    postfire: ee.Image,
    prefire: ee.Image,
    change: ee.Image,
    t_burn: float,
    t_rri_good: float,
    t_rri_low: float,
) -> ee.Image:
    rri = compute_rri(recent, postfire, prefire)
    dnbr = prefire.select("NBR").subtract(postfire.select("NBR"))
    burned = dnbr.gte(t_burn)
    ndre_good = recent.select("NDRE").gte(T_NDRE_GOOD)
    water = recent.select("NDWI").gt(T_WATER)
    valid = _recovery_valid_mask(recent, change)
    good = (rri.gte(t_rri_good).And(ndre_good).And(water.Not()).And(burned))
    failed = (rri.lt(t_rri_low).And(water.Not()).And(burned))
    cls = ee.Image.constant(2).rename("recovery_class")
    cls = cls.where(good, 1)
    cls = cls.where(failed, 3)
    cls = cls.where(valid.eq(1).And(water.Not()).And(burned.Not()), 4)
    cls = cls.where(valid.eq(0).Or(water), 0)
    return cls.clip(aoi).toUint8().rename("recovery_class")


def summarize_recovery_class(
    cls: ee.Image,
    parameter_name: str,
    parameter_value: float,
    t_burn: float,
    t_rri_good: float,
    t_rri_low: float,
    aoi: ee.Geometry,
    exp: dict,
    forest_mask_img,
) -> dict:
    if forest_mask_img is not None:
        cls = cls.updateMask(forest_mask_img)
    hist = cls.reduceRegion(
        reducer=ee.Reducer.frequencyHistogram(),
        geometry=aoi,
        scale=int(exp.get("scale", 20)),
        maxPixels=int(exp.get("max_pixels", 1e9)),
        bestEffort=True,
    ).getInfo()
    counts_raw = hist.get("recovery_class") or {}
    counts = {int(float(k)): int(v) for k, v in counts_raw.items()}
    c0 = counts.get(0, 0)
    c1 = counts.get(1, 0)
    c2 = counts.get(2, 0)
    c3 = counts.get(3, 0)
    c4 = counts.get(4, 0)
    total_px = c0 + c1 + c2 + c3 + c4
    burned_px = c1 + c2 + c3
    if total_px == 0: total_px = 1
    if burned_px == 0: burned_px = 1
    return {
        "parameter": parameter_name,
        "value": parameter_value,
        "burn_severity_threshold": t_burn,
        "rri_good_threshold": t_rri_good,
        "rri_low_threshold": t_rri_low,
        "class0_uncertain_px": c0,
        "class1_recovering_well_px": c1,
        "class2_recovering_weak_px": c2,
        "class3_failed_px": c3,
        "class4_outside_assessment_px": c4,
        "total_forest_mask_px": total_px,
        "burned_assessment_px": burned_px,
        "class0_uncertain_pct_total": 100.0 * c0 / total_px,
        "class1_recovering_well_pct_total": 100.0 * c1 / total_px,
        "class2_recovering_weak_pct_total": 100.0 * c2 / total_px,
        "class3_failed_pct_total": 100.0 * c3 / total_px,
        "class4_outside_assessment_pct_total": 100.0 * c4 / total_px,
        "recovering_well_pct_burned": 100.0 * c1 / burned_px,
        "recovering_weak_pct_burned": 100.0 * c2 / burned_px,
        "failed_pct_burned": 100.0 * c3 / burned_px,
    }
