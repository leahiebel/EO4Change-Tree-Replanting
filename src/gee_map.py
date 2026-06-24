"""
EO4Change Group 4 — Folium map builder with dynamic legend, Chart.js panel,
                     custom grouped layer panel, improved failed-cluster popups.
"""
import ee
import json
import base64
import folium
from pathlib import Path

from gee_config import (SRC_DIR, HTML_DIR, region_name, exp, aoi,
                         PREFIRE_YEARS, POSTFIRE_YEARS, RECENT_YEARS, CMP_MONTHS,
                         FM_ENABLED, FM_YEAR, FM_THRESHOLD,
                         T_RRI_GOOD, T_RRI_LOW, T_BSI_HIGH, T_WATER,
                         REF_DATE, REF_LABEL, ALL_BANDS,
                         SPEC_COLORS, BIOPH_COLORS, BIOPH_LABELS)
from gee_recovery import (RECOVERY_CLASSES, GUARD_CLASSES,
                            VIS_RECOVERY_CLASS, VIS_GUARDS)


# ── Legend helpers ────────────────────────────────────────────────────────────

def _safe_id(name: str) -> str:
    return "leg-" + "".join(c if c.isalnum() else "_" for c in name)


def _gradient_legend(title: str, palette: list, vmin: float, vmax: float,
                     unit: str = "", note: str = "") -> str:
    grad = ", ".join(palette)
    mid  = (vmin + vmax) / 2
    return (
        f'<div style="font-weight:600;margin-bottom:6px;font-size:13px;">{title}</div>'
        f'<div style="height:14px;background:linear-gradient(to right,{grad});'
        f'border:1px solid #aaa;border-radius:2px;margin-bottom:5px;"></div>'
        f'<div style="display:flex;justify-content:space-between;font-size:11px;margin-bottom:3px;">'
        f'  <span>{vmin:.2g}</span><span>{mid:.2g}</span><span>{vmax:.2g}</span>'
        f'</div>'
        + (f'<div style="font-size:11px;color:#555;">{unit}</div>' if unit else '')
        + (f'<div style="font-size:11px;color:#777;font-style:italic;margin-top:4px;">{note}</div>'
           if note else '')
    )


def _categorical_legend(title: str, classes: list, note: str = "") -> str:
    rows = "".join(
        f'<div style="display:flex;align-items:center;margin:3px 0;">'
        f'<span style="display:inline-block;width:16px;height:16px;background:{c};'
        f'border:1px solid #555;margin-right:7px;flex-shrink:0;"></span>'
        f'<span style="font-size:12px;">{idx} — {lbl}</span></div>'
        for idx, lbl, c in classes
    )
    return (
        f'<div style="font-weight:600;margin-bottom:6px;font-size:13px;">{title}</div>'
        + rows
        + (f'<div style="font-size:11px;color:#777;font-style:italic;margin-top:6px;">{note}</div>'
           if note else '')
    )


def _text_legend(title: str, note: str = "") -> str:
    return (
        f'<div style="font-weight:600;margin-bottom:5px;font-size:13px;">{title}</div>'
        f'<div style="font-size:11px;color:#777;font-style:italic;">'
        f'{note or "No radiometric scale."}</div>'
    )


def build_legends(lst_cfg, early_label, recent_label, postfire_label):
    """Build the _LEGENDS dict for the dynamic legend panel."""
    _change_bands_legends = {
        band: _gradient_legend(band, params["palette"], params["min"], params["max"],
                               note="Change = recent − post-fire composite")
        for band, params in {
            "NDVI_change": {"min": -0.3, "max": 0.3, "palette": ["#7d2222", "white", "#1f5e1f"]},
            "NDRE_change": {"min": -0.3, "max": 0.3, "palette": ["#7d2222", "white", "#1f5e1f"]},
            "NDMI_change": {"min": -0.3, "max": 0.3, "palette": ["#7d2222", "white", "#1f5e1f"]},
            "BSI_change":  {"min": -0.3, "max": 0.3, "palette": ["#1f5e1f", "white", "#7d2222"]},
            "NDWI_change": {"min": -0.3, "max": 0.3, "palette": ["saddlebrown", "white", "steelblue"]},
            "NBR_change":  {"min": -0.3, "max": 0.3, "palette": ["#7d2222", "white", "#1f5e1f"]},
        }.items()
    }

    legends = {
        "RGB":   _text_legend("True-colour composite (B4/B3/B2)"),
        "NDVI":  _gradient_legend("NDVI", ["saddlebrown", "khaki", "limegreen", "darkgreen"],
                                  -0.1, 0.85, note="Vegetation greenness — young forest ≈ 0.3–0.6"),
        "BSI":   _gradient_legend("BSI", ["darkgreen", "white", "orange", "saddlebrown"],
                                  -0.3, 0.4, note="High positive = exposed bare soil"),
        "NDWI":  _gradient_legend("NDWI", ["lightyellow", "powderblue", "steelblue", "navy"],
                                  -0.3, 0.5, note="Open water / flooding (McFeeters 1996)"),
        "NBR":   _gradient_legend("NBR", ["#7d2222", "#cc6600", "#ffe066", "#a8d666", "#1f5e1f"],
                                  -0.3, 0.7, note="Low = burned or stressed vegetation"),
        f"RRI (good≥{T_RRI_GOOD}, failed<{T_RRI_LOW})": _gradient_legend(
            "Relative Recovery Indicator (RRI)",
            ["#7d2222", "#ff6600", "#ffe066", "#a8d666", "#1f5e1f"],
            0, 1, note="(recent − postfire NDVI) / (prefire − postfire NDVI)"),
        "Recovery class": _categorical_legend(
            "Recovery class", RECOVERY_CLASSES,
            note="Spectral trajectory consistent with woody vegetation establishment"),
        "Uncertainty flags": _categorical_legend("Uncertainty flags", GUARD_CLASSES),
        early_label:   _text_legend(early_label,   f"Months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1} seasonal median"),
        postfire_label: _text_legend(postfire_label, f"Months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1} seasonal median"),
        recent_label:  _text_legend(recent_label,  f"Months {CMP_MONTHS[0]}–{CMP_MONTHS[1]-1} seasonal median"),
        **_change_bands_legends,
    }

    if FM_ENABLED:
        legends[f"Forest mask (DW {FM_YEAR}, p≥{FM_THRESHOLD})"] = _text_legend(
            f"Forest mask (DW {FM_YEAR}, p≥{FM_THRESHOLD})",
            "1 = was forest pre-fire (Dynamic World tree probability)"
        )

    if lst_cfg:
        legends.update({
            "LST median (Landsat)": _gradient_legend(
                "LST median (Landsat)",
                ["#313695", "#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#a50026"],
                15, 45, unit="°C",
                note=f"Target season {lst_cfg.get('target_year', '?')} — Landsat 8/9 C2 L2SP"),
            "LST anomaly (Landsat)": _gradient_legend(
                "LST anomaly (Landsat)",
                ["#313695", "#4575b4", "#e0f3f8", "#ffffff", "#fee090", "#a50026"],
                -5, 5, unit="°C", note="Target minus baseline median; blue = cooler than normal"),
            "LST stress class (Landsat)": _categorical_legend(
                "Thermal stress class",
                [(0, "No valid data", "#888888"), (1, "Near-normal", "#2ca25f"),
                 (2, "Moderate warm anomaly", "#feb24c"), (3, "Strong warm anomaly", "#de2d26")]),
        })

    return legends


def _cluster_reasoning(rows: list[dict], t_rri_low: float, t_bsi_high: float) -> str:
    """Generate reasoning bullets for why a cluster is class 3."""
    last = next((r for r in reversed(rows) if r.get("RRI") is not None), None)
    if last is None:
        return "Insufficient data to determine."
    bullets = []
    rri = last.get("RRI")
    bsi = last.get("BSI")
    if rri is not None and rri < t_rri_low:
        bullets.append(f"RRI = {rri:.3f} < {t_rri_low} (recovery ratio too low)")
    if bsi is not None and bsi > t_bsi_high:
        bullets.append(f"BSI = {bsi:.3f} > {t_bsi_high} (bare soil still elevated)")
    if not bullets:
        bullets.append("RRI/BSI values near thresholds — see time-series for context.")
    return "<br>".join(f"• {b}" for b in bullets)


def build_folium_map(
    aoi_bounds,
    raster,
    rri_img,
    recovery_class,
    guard_flags,
    forest_mask_img,
    prefire_composite,
    postfire_composite,
    recent_composite,
    change_img,
    VIS_CHANGE,
    lst_products,
    failed_diag_points,
    failed_diag_plot_paths,
    failed_diag_rows,
    records,
    lst_cfg,
) -> Path:
    """Build and save the Folium interactive map. Returns the output path."""

    _lons = [c[0] for c in aoi_bounds]
    _lats = [c[1] for c in aoi_bounds]
    _center = [(min(_lats) + max(_lats)) / 2, (min(_lons) + max(_lons)) / 2]

    fmap = folium.Map(location=_center, zoom_start=10, tiles="OpenStreetMap", control_scale=True)

    _visible_on_load: list = []
    # Maps display name → Folium JS variable name (e.g. "tile_abc123...")
    _layer_js_names: dict[str, str] = {}

    def _add_ee_layer(ee_image: ee.Image, params: dict, name: str, show: bool = False) -> None:
        if show:
            _visible_on_load.append(name)
        map_id = ee_image.getMapId(params)
        tile = folium.TileLayer(
            tiles=map_id["tile_fetcher"].url_format,
            attr="Google Earth Engine",
            name=name,
            overlay=True,
            control=True,
            show=show,
        )
        tile.add_to(fmap)
        _layer_js_names[name] = tile.get_name()

    # Composite label strings (used in legends and layer names)
    _RGB_PARAMS = {"min": 0, "max": 2800, "gamma": 1.4}
    early_label   = f"Pre-fire RGB ({PREFIRE_YEARS[0]} Apr–May)"
    postfire_label = f"Post-fire RGB ({POSTFIRE_YEARS[0]}–{POSTFIRE_YEARS[1]})"
    recent_label  = f"Recent RGB ({RECENT_YEARS[0]}–{RECENT_YEARS[1]})"
    _rri_layer_name = f"RRI (good≥{T_RRI_GOOD}, failed<{T_RRI_LOW})"

    VIS = {
        "RGB":  {"bands": ["B4", "B3", "B2"], "params": {"min": 0, "max": 2800, "gamma": 1.4}},
        "NDVI": {"bands": ["NDVI"], "params": {"min": -0.1, "max": 0.85,
                 "palette": ["saddlebrown", "khaki", "limegreen", "darkgreen"]}},
        "BSI":  {"bands": ["BSI"],  "params": {"min": -0.3, "max": 0.4,
                 "palette": ["darkgreen", "white", "orange", "saddlebrown"]}},
        "NDWI": {"bands": ["NDWI"], "params": {"min": -0.3, "max": 0.5,
                 "palette": ["lightyellow", "powderblue", "steelblue", "navy"]}},
        "NBR":  {"bands": ["NBR"],  "params": {"min": -0.3, "max": 0.7,
                 "palette": ["#7d2222", "#cc6600", "#ffe066", "#a8d666", "#1f5e1f"]}},
    }
    VIS_AVAILABLE = {k: v for k, v in VIS.items()
                     if all(b in ALL_BANDS + ["B4", "B3", "B2"] for b in v["bands"])}

    # Single-period summer raster layers
    for title, vc in VIS_AVAILABLE.items():
        _add_ee_layer(raster.select(vc["bands"]), vc["params"], title, show=(title == "RGB"))

    # RRI continuous layer
    _add_ee_layer(rri_img, {"min": 0, "max": 1, "palette": ["#7d2222", "#ff6600", "#ffe066", "#a8d666", "#1f5e1f"]},
                  _rri_layer_name, show=False)

    # Main recovery class map
    _add_ee_layer(recovery_class, VIS_RECOVERY_CLASS, "Recovery class", show=True)

    # Uncertainty / guard flags
    _add_ee_layer(guard_flags, VIS_GUARDS, "Uncertainty flags", show=False)

    # Forest mask
    _fm_key = f"Forest mask (DW {FM_YEAR}, p≥{FM_THRESHOLD})"
    if forest_mask_img is not None:
        _add_ee_layer(forest_mask_img, {"min": 0, "max": 1, "palette": ["#dddddd", "#1a6b1a"]},
                      _fm_key, show=False)

    # Pre-fire, post-fire, recent RGB composites
    _add_ee_layer(prefire_composite.select(["B4", "B3", "B2"]),  _RGB_PARAMS, early_label,   show=False)
    _add_ee_layer(postfire_composite.select(["B4", "B3", "B2"]), _RGB_PARAMS, postfire_label, show=False)
    _add_ee_layer(recent_composite.select(["B4", "B3", "B2"]),   _RGB_PARAMS, recent_label,  show=False)

    # Change bands
    for _band, _params in VIS_CHANGE.items():
        _add_ee_layer(change_img.select(_band), _params, _band, show=False)

    # LST layers
    if lst_products:
        _LST_PALETTE  = ["#313695", "#4575b4", "#abd9e9", "#ffffbf", "#fdae61", "#a50026"]
        _ANOM_PALETTE = ["#313695", "#4575b4", "#e0f3f8", "#ffffff", "#fee090", "#a50026"]
        _STRESS_PALETTE = ["#888888", "#2ca25f", "#feb24c", "#de2d26"]
        _add_ee_layer(lst_products["lst_median"],
                      {"min": 15, "max": 45, "palette": _LST_PALETTE},
                      "LST median (Landsat)", show=False)
        _add_ee_layer(lst_products["lst_anomaly"],
                      {"min": -5, "max": 5, "palette": _ANOM_PALETTE},
                      "LST anomaly (Landsat)", show=False)
        _add_ee_layer(lst_products["lst_stress"],
                      {"min": 0, "max": 3, "palette": _STRESS_PALETTE},
                      "LST stress class (Landsat)", show=False)

    # AOI outline — capture JS name for custom panel
    _aoi_geojson = folium.GeoJson(
        aoi.getInfo(),
        name="AOI",
        style_function=lambda _f: {"color": "red", "weight": 2, "fillOpacity": 0},
    )
    _aoi_geojson.add_to(fmap)
    _layer_js_names["AOI"] = _aoi_geojson.get_name()
    _visible_on_load.append("AOI")

    # ── Build cluster_ts_json for CLUSTER_DATA JS variable ───────────────────
    cluster_ts_json: dict = {}
    for row in failed_diag_rows:
        eid = row["example_id"]
        if eid not in cluster_ts_json:
            cluster_ts_json[eid] = []
        cluster_ts_json[eid].append({
            k: row.get(k)
            for k in ["date", "NDVI", "RRI", "NDRE", "BSI", "NDMI", "NDWI", "NBR"]
        })

    # ── Improved failed cluster markers ──────────────────────────────────────
    if failed_diag_points:
        failed_group = folium.FeatureGroup(name="Failed-pixel diagnostics", show=True)

        for ex in failed_diag_points:
            ex_id     = ex["id"]
            plot_path = failed_diag_plot_paths.get(ex_id)

            ex_rows = [r for r in failed_diag_rows if r["example_id"] == ex_id]
            rank    = next((r.get("cluster_rank", "?") for r in ex_rows), "?")

            reasoning_html = _cluster_reasoning(ex_rows, T_RRI_LOW, T_BSI_HIGH)

            png_img_tag = ""
            if plot_path is not None and Path(plot_path).exists():
                with open(plot_path, "rb") as f:
                    img64 = base64.b64encode(f.read()).decode("utf-8")
                png_img_tag = (
                    f'<img src="data:image/png;base64,{img64}" '
                    f'width="660" style="margin-top:8px;display:block;">'
                )

            _has_ts = ex_id in cluster_ts_json
            chart_btn = (
                f'<div style="margin-bottom:8px;">'
                f'<button onclick="showClusterInPanel(\'{ex_id}\')"'
                f' style="background:#1f5e1f;color:white;border:none;border-radius:3px;'
                f'padding:4px 10px;cursor:pointer;font-size:11px;">'
                f'Open interactive time-series →</button></div>'
            ) if _has_ts else ""

            popup_html = f"""
<div style="width:680px;">
  <h4 style="margin:4px 0 4px 0;">{ex_id} — Cluster rank #{rank}</h4>
  <div style="font-size:11px;color:#555;margin-bottom:8px;">
    lat {ex['lat']:.5f}&nbsp;&nbsp;lon {ex['lon']:.5f}
  </div>
  <div style="background:#fff3cd;padding:6px 8px;border-radius:3px;
              font-size:11px;margin-bottom:8px;line-height:1.6;">
    <b>Why class 3 — Not recovering:</b><br>{reasoning_html}
  </div>
  {chart_btn}
  {png_img_tag}
</div>
"""
            folium.CircleMarker(
                location=[ex["lat"], ex["lon"]],
                radius=7,
                color="#ff0000",
                fill=True,
                fill_color="#ff0000",
                fill_opacity=0.9,
                popup=folium.Popup(popup_html, max_width=720),
                tooltip=f"{ex_id} — failed recovery (click for details)",
            ).add_to(failed_group)

        failed_group.add_to(fmap)
        _layer_js_names["Failed-pixel diagnostics"] = failed_group.get_name()
        _visible_on_load.append("Failed-pixel diagnostics")

    # ── Build _LEGENDS dict ───────────────────────────────────────────────────
    _LEGENDS = build_legends(lst_cfg, early_label, recent_label, postfire_label)

    # ── Dynamic legend (bottom-right, slightly larger) ────────────────────────
    _inner_divs = "\n".join(
        f'<div id="{_safe_id(n)}" '
        f'style="display:{"block" if n in _visible_on_load else "none"};">'
        f'{html}</div>'
        for n, html in _LEGENDS.items()
    )
    _any_visible_init = "none" if _visible_on_load else "block"
    _legend_outer = (
        '<div id="dyn-legend" style="'
        'position:fixed;bottom:30px;right:30px;z-index:9999;'
        'background:rgba(255,255,255,0.97);padding:13px 15px;'
        'border:1px solid #888;border-radius:5px;'
        'font:13px/1.4 system-ui,sans-serif;max-width:300px;'
        'box-shadow:0 2px 8px rgba(0,0,0,0.22);">'
        + _inner_divs
        + f'<div id="dyn-legend-empty" style="display:{_any_visible_init};font-size:12px;color:#888;">'
          f'Toggle a layer to see its legend.</div>'
        + '</div>'
    )
    fmap.get_root().html.add_child(folium.Element(_legend_outer))

    _id_map_js = json.dumps({n: _safe_id(n) for n in _LEGENDS})
    _map_var   = fmap.get_name()
    _dyn_script = f"""
(function poll() {{
  var m = window['{_map_var}'];
  if (!m) {{ setTimeout(poll, 50); return; }}
  var ID_MAP = {_id_map_js};
  function refresh() {{
    var any = Object.values(ID_MAP).some(function(id) {{
      var el = document.getElementById(id);
      return el && el.style.display !== 'none';
    }});
    var emp = document.getElementById('dyn-legend-empty');
    if (emp) emp.style.display = any ? 'none' : 'block';
  }}
  function setLayer(name, visible) {{
    var id = ID_MAP[name];
    if (!id) return;
    var el = document.getElementById(id);
    if (el) el.style.display = visible ? 'block' : 'none';
    refresh();
  }}
  m.on('overlayadd',    function(e) {{ setLayer(e.name, true);  }});
  m.on('overlayremove', function(e) {{ setLayer(e.name, false); }});
  refresh();
}})();
"""
    fmap.get_root().script.add_child(folium.Element(_dyn_script))

    # ── Chart.js side panel ───────────────────────────────────────────────────
    fmap.get_root().header.add_child(folium.Element(
        '<script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>\n'
        '<script src="https://cdn.jsdelivr.net/npm/chartjs-adapter-date-fns@3.0.0/dist/chartjs-adapter-date-fns.bundle.min.js"></script>\n'
        '<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@3.0.1/dist/chartjs-plugin-annotation.min.js"></script>'
    ))

    _ts_records_clean = [
        {k: r.get(k) for k in ["date", "label"] + ALL_BANDS + ["RRI"]}
        for r in records
    ]
    _ts_records_json = json.dumps(_ts_records_clean)

    _ts_layer_vars = {
        "RGB":                  ["NDVI"],
        "NDVI":                 ["NDVI"],
        "BSI":                  ["BSI"],
        "NDWI":                 ["NDWI"],
        "NBR":                  ["NBR"],
        _rri_layer_name:        ["RRI"],
        "Recovery class":       ["NDVI", "BSI", "NDRE"],
        "Uncertainty flags":    ["NDVI", "BSI"],
        early_label:            ["NDVI"],
        postfire_label:         ["NDVI"],
        recent_label:           ["NDVI"],
        "NDVI_change":          ["NDVI"],
        "NDRE_change":          ["NDRE"],
        "NDMI_change":          ["NDMI"],
        "BSI_change":           ["BSI"],
        "NDWI_change":          ["NDWI"],
        "NBR_change":           ["NBR"],
    }
    _ts_layer_vars_json = json.dumps(_ts_layer_vars)
    _ts_colors_json     = json.dumps({**SPEC_COLORS, **BIOPH_COLORS, "RRI": "#9b59b6"})
    _ts_ref_date        = json.dumps(REF_DATE  or "")
    _ts_ref_label       = json.dumps(REF_LABEL or "")
    _ts_init_json       = json.dumps(_visible_on_load)
    _cluster_data_json  = json.dumps(cluster_ts_json)
    _rri_low_js         = json.dumps(T_RRI_LOW)
    _bsi_high_js        = json.dumps(T_BSI_HIGH)

    _ts_panel_html = (
        '<div id="ts-panel" style="'
        'position:fixed;left:50px;top:60px;z-index:9998;'
        'background:rgba(255,255,255,0.97);'
        'padding:12px 14px;border:1px solid #888;border-radius:4px;'
        'box-shadow:0 2px 8px rgba(0,0,0,0.25);width:340px;'
        'font:11px/1.4 system-ui,sans-serif;display:none;">'
        '<div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:6px;">'
        '<span id="ts-title" style="font-weight:600;font-size:12px;max-width:280px;'
        'overflow:hidden;text-overflow:ellipsis;white-space:nowrap;"></span>'
        '<button id="ts-close" style="background:none;border:none;cursor:pointer;'
        'font-size:18px;line-height:1;color:#555;flex-shrink:0;">×</button>'
        '</div>'
        '<div style="position:relative;height:200px;">'
        '<canvas id="ts-chart"></canvas>'
        '</div>'
        '</div>'
    )
    fmap.get_root().html.add_child(folium.Element(_ts_panel_html))

    _ts_script = f"""
(function waitForChart() {{
  if (typeof Chart === 'undefined') {{ setTimeout(waitForChart, 50); return; }}

  var TS_DATA         = {_ts_records_json};
  var LAYER_VARS      = {_ts_layer_vars_json};
  var COLORS          = {_ts_colors_json};
  var REF_DATE        = {_ts_ref_date};
  var REF_LABEL       = {_ts_ref_label};
  var INIT_LAYERS     = {_ts_init_json};
  var CLUSTER_DATA    = {_cluster_data_json};
  var RRI_LOW_THRESHOLD = {_rri_low_js};
  var BSI_HIGH_THRESHOLD = {_bsi_high_js};

  var panel   = document.getElementById('ts-panel');
  var tsChart = null;

  document.getElementById('ts-close').addEventListener('click', function() {{
    panel.style.display = 'none';
  }});

  function showPanel(layerName) {{
    var vars = LAYER_VARS[layerName];
    if (!vars || vars.length === 0) return;
    var pts_check = TS_DATA.filter(function(r) {{
      return vars.some(function(v) {{ return r[v] !== null && r[v] !== undefined; }});
    }});
    if (pts_check.length === 0) return;
    document.getElementById('ts-title').textContent = layerName;
    panel.style.display = 'block';
    var datasets = vars.map(function(v) {{
      var pts = TS_DATA
        .filter(function(r) {{ return r[v] !== null && r[v] !== undefined; }})
        .map(function(r)   {{ return {{ x: r.date, y: r[v] }}; }});
      return {{
        label: v,
        data: pts,
        borderColor: COLORS[v] || '#4a9',
        backgroundColor: (COLORS[v] || '#4a9') + '22',
        borderWidth: 1.5,
        pointRadius: 2.5,
        pointHoverRadius: 4,
        tension: 0.2,
        fill: false,
      }};
    }});
    var annotations = {{}};
    if (REF_DATE) {{
      annotations['ref'] = {{
        type: 'line', xMin: REF_DATE, xMax: REF_DATE,
        borderColor: 'rgba(210,40,40,0.8)', borderWidth: 1.5, borderDash: [6, 4],
        label: {{
          display: true, content: REF_LABEL, position: 'start',
          backgroundColor: 'rgba(210,40,40,0.1)', color: '#c00', font: {{ size: 9 }}, padding: 3,
        }},
      }};
    }}
    if (tsChart) tsChart.destroy();
    var ctx = document.getElementById('ts-chart').getContext('2d');
    tsChart = new Chart(ctx, {{
      type: 'line',
      data: {{ datasets: datasets }},
      options: {{
        responsive: true, maintainAspectRatio: false,
        interaction: {{ mode: 'index', intersect: false }},
        scales: {{
          x: {{
            type: 'time', time: {{ unit: 'year', tooltipFormat: 'yyyy-MM-dd' }},
            ticks: {{ maxTicksLimit: 6, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
          y: {{
            ticks: {{ maxTicksLimit: 5, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
        }},
        plugins: {{
          legend: {{ display: vars.length > 1, labels: {{ font: {{ size: 9 }}, boxWidth: 12 }} }},
          tooltip: {{
            callbacks: {{
              label: function(ctx) {{
                var v = ctx.parsed.y;
                return ctx.dataset.label + ': ' + (v !== null ? v.toFixed(3) : 'N/A');
              }},
            }},
          }},
          annotation: {{ annotations: annotations }},
        }},
      }},
    }});
  }}

  window.showClusterInPanel = function(clusterId) {{
    var clusterData = CLUSTER_DATA[clusterId];
    if (!clusterData) return;
    var vars = ['NDVI', 'RRI', 'NDRE', 'BSI'];
    var datasets = vars.map(function(v) {{
      var pts = clusterData
        .filter(function(r) {{ return r[v] != null; }})
        .map(function(r)   {{ return {{ x: r.date, y: r[v] }}; }});
      return {{
        label: v,
        data: pts,
        borderColor: COLORS[v] || '#888',
        backgroundColor: (COLORS[v] || '#888') + '22',
        borderWidth: 1.5, pointRadius: 2.5, tension: 0.2, fill: false,
      }};
    }});
    var annotations = {{}};
    annotations['rri_low'] = {{
      type: 'line', yMin: RRI_LOW_THRESHOLD, yMax: RRI_LOW_THRESHOLD,
      borderColor: '#cc0000', borderWidth: 1, borderDash: [4, 3],
      label: {{
        display: true, content: 'RRI failed threshold', position: 'end',
        backgroundColor: 'rgba(200,0,0,0.08)', color: '#cc0000',
        font: {{ size: 8 }}, padding: 2,
      }},
    }};
    annotations['bsi_high'] = {{
      type: 'line', yMin: BSI_HIGH_THRESHOLD, yMax: BSI_HIGH_THRESHOLD,
      borderColor: '#c47a1e', borderWidth: 1, borderDash: [4, 3],
      label: {{
        display: true, content: 'BSI high threshold', position: 'end',
        backgroundColor: 'rgba(196,122,30,0.08)', color: '#c47a1e',
        font: {{ size: 8 }}, padding: 2,
      }},
    }};
    if (REF_DATE) {{
      annotations['ref'] = {{
        type: 'line', xMin: REF_DATE, xMax: REF_DATE,
        borderColor: 'rgba(210,40,40,0.8)', borderWidth: 1.5, borderDash: [6, 4],
        label: {{
          display: true, content: REF_LABEL, position: 'start',
          backgroundColor: 'rgba(210,40,40,0.1)', color: '#c00', font: {{ size: 9 }}, padding: 3,
        }},
      }};
    }}
    document.getElementById('ts-title').textContent = clusterId + ' — failed cluster';
    panel.style.display = 'block';
    if (tsChart) tsChart.destroy();
    var ctx = document.getElementById('ts-chart').getContext('2d');
    tsChart = new Chart(ctx, {{
      type: 'line',
      data: {{ datasets: datasets }},
      options: {{
        responsive: true, maintainAspectRatio: false,
        interaction: {{ mode: 'index', intersect: false }},
        scales: {{
          x: {{
            type: 'time', time: {{ unit: 'year', tooltipFormat: 'yyyy-MM-dd' }},
            ticks: {{ maxTicksLimit: 6, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
          y: {{
            ticks: {{ maxTicksLimit: 5, font: {{ size: 9 }} }},
            grid: {{ color: 'rgba(0,0,0,0.06)' }},
          }},
        }},
        plugins: {{
          legend: {{ display: true, labels: {{ font: {{ size: 9 }}, boxWidth: 12 }} }},
          annotation: {{ annotations: annotations }},
        }},
      }},
    }});
  }};

  (function pollMap() {{
    var m = window['{_map_var}'];
    if (!m) {{ setTimeout(pollMap, 50); return; }}
    m.on('overlayadd', function(e) {{ showPanel(e.name); }});
    m.on('overlayremove', function(e) {{
      var title = document.getElementById('ts-title');
      if (title && title.textContent === e.name) {{ panel.style.display = 'none'; }}
    }});
    for (var i = 0; i < INIT_LAYERS.length; i++) {{
      if (LAYER_VARS[INIT_LAYERS[i]]) {{
        showPanel(INIT_LAYERS[i]);
        break;
      }}
    }}
  }})();
}})();
"""
    fmap.get_root().script.add_child(folium.Element(_ts_script))

    # ── Custom grouped layer panel (replaces folium.LayerControl) ────────────

    # Main product: 3 layers with opacity sliders
    _main_product = [
        ("Recovery class",    "Recovery class",  True,  "#1f5e1f"),
        ("RRI",               _rri_layer_name,   False, "#9b59b6"),
        ("Uncertainty flags", "Uncertainty flags", False, "#00b894"),
    ]

    # Other-layer groups (built dynamically from what was actually added)
    _other_groups: list[tuple[str, list[tuple[str, bool]]]] = []

    _idx_layers = [(n, n in _visible_on_load)
                   for n in ["RGB", "NDVI", "BSI", "NDWI", "NBR"]
                   if n in _layer_js_names]
    if _idx_layers:
        _other_groups.append(("Spectral indices", _idx_layers))

    if FM_ENABLED and _fm_key in _layer_js_names:
        _other_groups.append(("Forest mask", [(_fm_key, False)]))

    _comp_layers = [(n, False) for n in [early_label, postfire_label, recent_label]
                    if n in _layer_js_names]
    if _comp_layers:
        _other_groups.append(("Composites", _comp_layers))

    _chg_layers = [(n, False) for n in
                   ["NDVI_change", "NDRE_change", "NDMI_change",
                    "BSI_change", "NDWI_change", "NBR_change"]
                   if n in _layer_js_names]
    if _chg_layers:
        _other_groups.append(("Change bands", _chg_layers))

    if lst_products:
        _lst_layers = [(n, False) for n in
                       ["LST median (Landsat)", "LST anomaly (Landsat)", "LST stress class (Landsat)"]
                       if n in _layer_js_names]
        if _lst_layers:
            _other_groups.append(("Land Surface Temp.", _lst_layers))

    _diag_layers = [(n, n in _visible_on_load)
                    for n in ["AOI", "Failed-pixel diagnostics"]
                    if n in _layer_js_names]
    if _diag_layers:
        _other_groups.append(("Diagnostics", _diag_layers))

    # Build the "Recovery Products" section HTML
    def _pct_id(label: str) -> str:
        return "pct-" + "".join(c if c.isalnum() else "_" for c in label)

    def _chk_id(label: str) -> str:
        return "chk-" + "".join(c if c.isalnum() else "_" for c in label)

    _main_rows_html = ""
    for short_label, full_name, initial_show, accent in _main_product:
        checked = "checked" if initial_show else ""
        _main_rows_html += f"""
<div style="margin-bottom:10px;">
  <label style="display:flex;align-items:center;gap:7px;cursor:pointer;font-size:13px;font-weight:500;">
    <input type="checkbox" id="{_chk_id(full_name)}" data-lname="{full_name}"
           onchange="lpToggle(this)" {checked}
           style="width:15px;height:15px;cursor:pointer;accent-color:{accent};">
    {short_label}
  </label>
  <div style="display:flex;align-items:center;gap:6px;margin-top:4px;padding-left:22px;">
    <input type="range" min="0" max="100" value="100" data-lname="{full_name}"
           oninput="lpOpacity(this)"
           style="flex:1;height:4px;cursor:pointer;accent-color:{accent};">
    <span id="{_pct_id(full_name)}" style="font-size:10px;color:#777;width:30px;text-align:right;">100%</span>
  </div>
</div>"""

    # Build the "Other layers" groups HTML
    _other_html = ""
    for grp_title, grp_layers in _other_groups:
        _other_html += (
            f'<div style="margin-bottom:8px;">'
            f'<div style="font-size:10px;font-weight:700;color:#888;text-transform:uppercase;'
            f'letter-spacing:0.05em;margin-bottom:5px;padding-bottom:3px;border-bottom:1px solid #eee;">'
            f'{grp_title}</div>'
        )
        for name, init_show in grp_layers:
            checked = "checked" if init_show else ""
            _other_html += (
                f'<label style="display:flex;align-items:center;gap:6px;cursor:pointer;'
                f'font-size:12px;margin-bottom:4px;">'
                f'<input type="checkbox" data-lname="{name}" onchange="lpToggle(this)" {checked}'
                f' style="width:13px;height:13px;cursor:pointer;">'
                f'{name}</label>'
            )
        _other_html += '</div>'

    _panel_html = f"""
<div id="lp-panel" style="
  position:fixed; top:80px; right:10px; z-index:9999;
  background:rgba(255,255,255,0.97);
  border:1px solid #bbb; border-radius:6px;
  box-shadow:0 3px 10px rgba(0,0,0,0.2);
  font:13px/1.4 system-ui,sans-serif;
  width:250px; user-select:none; overflow:hidden;">

  <!-- Header -->
  <div style="background:#1f5e1f;color:white;padding:9px 13px;
              font-weight:700;font-size:13px;letter-spacing:0.02em;">
    Layers
  </div>

  <!-- Recovery products section -->
  <div style="padding:12px 13px;border-bottom:1px solid #e0e0e0;">
    <div style="font-size:10px;font-weight:700;color:#1f5e1f;text-transform:uppercase;
                letter-spacing:0.07em;margin-bottom:10px;">
      ● Recovery products
    </div>
    {_main_rows_html}
  </div>

  <!-- Other layers (collapsible) -->
  <div>
    <div id="lp-other-hdr" onclick="lpToggleOther()"
         style="padding:9px 13px;cursor:pointer;display:flex;
                justify-content:space-between;align-items:center;
                background:#f7f7f7;border-bottom:1px solid #e0e0e0;">
      <span style="font-size:10px;font-weight:700;color:#666;
                   text-transform:uppercase;letter-spacing:0.07em;">Other layers</span>
      <span id="lp-arrow" style="color:#888;font-size:11px;">▶</span>
    </div>
    <div id="lp-other-body" style="display:none;padding:10px 13px;
                                   max-height:340px;overflow-y:auto;">
      {_other_html}
    </div>
  </div>
</div>
"""
    fmap.get_root().html.add_child(folium.Element(_panel_html))

    # JS variable name → actual JS object map (passed to JS as a JSON dict)
    _js_name_map_json = json.dumps(_layer_js_names)

    _panel_script = f"""
(function waitForMap() {{
  var m = window['{_map_var}'];
  if (!m) {{ setTimeout(waitForMap, 50); return; }}

  var LP_JS_MAP = {_js_name_map_json};

  function lpGetLayer(name) {{
    var jsVar = LP_JS_MAP[name];
    return jsVar ? window[jsVar] : null;
  }}

  window.lpToggle = function(checkbox) {{
    var name  = checkbox.getAttribute('data-lname');
    var layer = lpGetLayer(name);
    if (!layer) return;
    if (checkbox.checked) {{
      layer.addTo(m);
      m.fire('overlayadd',    {{name: name, layer: layer}});
    }} else {{
      m.removeLayer(layer);
      m.fire('overlayremove', {{name: name, layer: layer}});
    }}
  }};

  window.lpOpacity = function(slider) {{
    var name  = slider.getAttribute('data-lname');
    var layer = lpGetLayer(name);
    if (layer && layer.setOpacity) layer.setOpacity(slider.value / 100);
    var pctEl = document.getElementById('pct-' + name.replace(/[^a-zA-Z0-9]/g, '_'));
    if (pctEl) pctEl.textContent = slider.value + '%';
  }};

  window.lpToggleOther = function() {{
    var body  = document.getElementById('lp-other-body');
    var arrow = document.getElementById('lp-arrow');
    if (body.style.display === 'none') {{
      body.style.display = 'block';
      arrow.textContent  = '▼';
    }} else {{
      body.style.display = 'none';
      arrow.textContent  = '▶';
    }}
  }};
}})();
"""
    fmap.get_root().script.add_child(folium.Element(_panel_script))

    out_map = HTML_DIR / f"map_{region_name}.html"
    fmap.save(str(out_map))
    return out_map
