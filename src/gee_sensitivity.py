"""
EO4Change Group 4 — OAT threshold sensitivity analysis.
"""
import pandas as pd
import matplotlib.pyplot as plt
from gee_config import (SRC_DIR, FIGURES_DIR, CSV_DIR,
                         region_name, exp, aoi,
                         T_BURN, T_RRI_GOOD, T_RRI_LOW,
                         SENS_BURN_VALUES, SENS_RRI_GOOD_VALUES, SENS_RRI_LOW_VALUES)
from gee_recovery import build_recovery_class_sensitivity, summarize_recovery_class


def run_threshold_sensitivity_analysis(
    recent_composite,
    postfire_composite,
    prefire_composite,
    change_img,
    forest_mask_img,
) -> None:
    rows = []
    print("Running threshold sensitivity analysis...")

    baseline_cls = build_recovery_class_sensitivity(
        recent=recent_composite, postfire=postfire_composite,
        prefire=prefire_composite, change=change_img,
        t_burn=T_BURN, t_rri_good=T_RRI_GOOD, t_rri_low=T_RRI_LOW,
    )
    rows.append(summarize_recovery_class(
        cls=baseline_cls, parameter_name="baseline", parameter_value=0.0,
        t_burn=T_BURN, t_rri_good=T_RRI_GOOD, t_rri_low=T_RRI_LOW,
        aoi=aoi, exp=exp, forest_mask_img=forest_mask_img,
    ))

    for value in SENS_BURN_VALUES:
        value = float(value)
        print(f"  testing burn_severity_threshold = {value}")
        cls = build_recovery_class_sensitivity(
            recent=recent_composite, postfire=postfire_composite,
            prefire=prefire_composite, change=change_img,
            t_burn=value, t_rri_good=T_RRI_GOOD, t_rri_low=T_RRI_LOW,
        )
        rows.append(summarize_recovery_class(
            cls=cls, parameter_name="burn_severity_threshold", parameter_value=value,
            t_burn=value, t_rri_good=T_RRI_GOOD, t_rri_low=T_RRI_LOW,
            aoi=aoi, exp=exp, forest_mask_img=forest_mask_img,
        ))

    for value in SENS_RRI_GOOD_VALUES:
        value = float(value)
        print(f"  testing rri_good_threshold = {value}")
        cls = build_recovery_class_sensitivity(
            recent=recent_composite, postfire=postfire_composite,
            prefire=prefire_composite, change=change_img,
            t_burn=T_BURN, t_rri_good=value, t_rri_low=T_RRI_LOW,
        )
        rows.append(summarize_recovery_class(
            cls=cls, parameter_name="rri_good_threshold", parameter_value=value,
            t_burn=T_BURN, t_rri_good=value, t_rri_low=T_RRI_LOW,
            aoi=aoi, exp=exp, forest_mask_img=forest_mask_img,
        ))

    for value in SENS_RRI_LOW_VALUES:
        value = float(value)
        print(f"  testing rri_low_threshold = {value}")
        cls = build_recovery_class_sensitivity(
            recent=recent_composite, postfire=postfire_composite,
            prefire=prefire_composite, change=change_img,
            t_burn=T_BURN, t_rri_good=T_RRI_GOOD, t_rri_low=value,
        )
        rows.append(summarize_recovery_class(
            cls=cls, parameter_name="rri_low_threshold", parameter_value=value,
            t_burn=T_BURN, t_rri_good=T_RRI_GOOD, t_rri_low=value,
            aoi=aoi, exp=exp, forest_mask_img=forest_mask_img,
        ))

    df = pd.DataFrame(rows)
    baseline = df[df["parameter"] == "baseline"].iloc[0]
    for col in ["recovering_well_pct_burned", "recovering_weak_pct_burned",
                "failed_pct_burned", "class4_outside_assessment_pct_total", "burned_assessment_px"]:
        df[f"delta_{col}"] = df[col] - baseline[col]

    out_csv = CSV_DIR / f"sensitivity_thresholds_{region_name}.csv"
    df.to_csv(out_csv, index=False)
    print(f"Sensitivity CSV → {out_csv}")

    for parameter_name in ["burn_severity_threshold", "rri_good_threshold", "rri_low_threshold"]:
        sub = df[df["parameter"] == parameter_name].copy()
        if sub.empty:
            continue
        sub = sub.sort_values("value")
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(sub["value"], sub["recovering_well_pct_burned"], marker="o",
                label="Recovering well (% of burned assessment)")
        ax.plot(sub["value"], sub["recovering_weak_pct_burned"], marker="o",
                label="Recovering weak (% of burned assessment)")
        ax.plot(sub["value"], sub["failed_pct_burned"], marker="o",
                label="Failed (% of burned assessment)")
        if parameter_name == "burn_severity_threshold":
            ax.plot(sub["value"], sub["class4_outside_assessment_pct_total"], marker="o",
                    linestyle="--", label="Outside assessment (% of total forest mask)")
        ax.axvline(
            {"burn_severity_threshold": T_BURN, "rri_good_threshold": T_RRI_GOOD,
             "rri_low_threshold": T_RRI_LOW}[parameter_name],
            linestyle=":", linewidth=1.5, label="Baseline value",
        )
        ax.set_xlabel(parameter_name)
        ax.set_ylabel("Pixel percentage")
        ax.set_title(f"Sensitivity analysis: {parameter_name}")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
        plt.tight_layout()
        out_png = FIGURES_DIR / f"sensitivity_{parameter_name}_{region_name}.png"
        fig.savefig(out_png, dpi=150)
        plt.close(fig)
        print(f"Sensitivity plot → {out_png}")

    print("\nSensitivity summary:")
    cols_to_print = ["parameter", "value", "recovering_well_pct_burned",
                     "recovering_weak_pct_burned", "failed_pct_burned",
                     "class4_outside_assessment_pct_total"]
    print(df[cols_to_print].to_string(index=False))
