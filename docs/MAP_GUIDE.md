# Reading the EO4Change interactive map

A field guide to `src/map_<region>.html` — what each layer is, how it was
computed, what its colors mean, and how to interpret what you see.

The goal of this document is to let a reader who has *never seen this
project before* understand, in one sitting, exactly what they are looking
at on the map and why we drew it that way.

---

## 1. Big picture — what are we trying to do?

We want to know, **per pixel**, whether a piece of land that was
replanted (or that we expect to be revegetating) is on a healthy
trajectory, stagnating, or struggling.

To do that we compare two satellite snapshots of the same place:

| Period | Role | Default for Portugal |
|---|---|---|
| **Early** composite | "before" — the reference state | summers 2017–2018 |
| **Recent** composite | "after" — the current state | summers 2022–2025 |

Each composite is a **median image** built from many Sentinel-2 scenes
inside a fixed seasonal window (May–August by default). Taking the
median over multiple years and months suppresses individual cloudy
days, water-vapour artefacts, and one-off pixel noise — what remains is
a stable "typical summer" image for that period.

Then we do three things, in this order:

1. **Single-period view** — show the recent summer image as natural
   color, NDVI, NBR, etc. so you can see *what is there now*.
2. **Change layers** — subtract early from recent for each index, so
   you can see *what got better or worse*.
3. **Establishment status map** — collapse those change layers into
   four discrete classes (`Good / Moderate / Weak / Uncertain`) using
   explicit rules.

Each of those three sits in the layer panel (top right of the map) and
can be toggled independently.

---

## 2. Where the pixels come from

### 2.1 Satellite

**Sentinel-2 L2A "Harmonized"** (`COPERNICUS/S2_SR_HARMONIZED` in
Earth Engine). Surface reflectance, 10–20 m per pixel depending on
band, revisit ≈ 5 days at the equator (faster in Europe with two
satellites). "Harmonized" means ESA's processing change in
January 2022 has been compensated — pre-2022 and post-2022 scenes are
on the same radiometric scale, which matters a lot for change
detection.

### 2.2 Cloud masking

We use **Cloud Score Plus** (`GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED`).
Each pixel gets a "clear-sky" score between 0 (definitely cloud /
shadow / haze) and 1 (clean). We keep pixels with score ≥ **0.70**
(configurable via `options.cs_plus_threshold`).

Cloud Score Plus is a learned model trained on hand-labelled scenes,
which makes it noticeably more accurate than the standard SCL mask —
especially around thin cirrus and shadow edges, which are exactly the
pixels that confuse vegetation indices.

### 2.3 Compositing

For each composite period we:
1. Filter S2 by date window and AOI.
2. Apply the Cloud Score Plus mask.
3. Compute the spectral indices on the *masked* surface reflectance.
4. Take the **per-pixel median across all surviving scenes**.

Median (rather than mean) is robust to the cloud edges that occasionally
sneak past the mask: a single anomalous value cannot drag the result.

---

## 3. The spectral indices

All Sentinel-2 indices below use band numbers from the satellite's
official band table. Wavelengths are in nanometres.

| Band | Wavelength | Native res | What it sees |
|---|---|---|---|
| B2  | 490 nm  | 10 m | Blue |
| B3  | 560 nm  | 10 m | Green |
| B4  | 665 nm  | 10 m | Red — strongly absorbed by chlorophyll |
| B5  | 705 nm  | 20 m | Red-edge 1 |
| B8  | 842 nm  | 10 m | Near-infrared (NIR) — strongly reflected by healthy leaves |
| B8A | 865 nm  | 20 m | Narrow NIR |
| B11 | 1610 nm | 20 m | Shortwave infrared 1 (SWIR1) — sensitive to leaf and soil water |
| B12 | 2190 nm | 20 m | SWIR2 — sensitive to dry biomass, burn scars |

The physical logic behind every vegetation index is the same: healthy
leaves *reflect a lot of NIR and absorb red*, while bare soil reflects
roughly equally in red and NIR, and water absorbs almost all NIR.
Different bands isolate different stress signals (water content, lignin,
chlorophyll concentration, etc.).

### 3.1 NDVI — Normalised Difference Vegetation Index

$$
\mathrm{NDVI} = \frac{B8 - B4}{B8 + B4}
$$

- **What it measures:** density and greenness of photosynthetically
  active vegetation. The single most used vegetation index in remote
  sensing.
- **Range:** −1 to +1. In practice: water ≈ −0.1 to 0.0, bare soil ≈
  0.1–0.2, sparse grass ≈ 0.3, healthy crop/forest ≈ 0.6–0.9.
- **Color ramp on the map:** brown → yellow → green (low → high).
  Bright green = vigorous canopy.
- **Failure modes:** saturates over dense canopy (an old-growth forest
  and a 10-year-old plantation can both read ≈ 0.85), and is fooled by
  bright understorey grass on a thinly stocked plantation.

### 3.2 NDRE — Normalised Difference Red-Edge

$$
\mathrm{NDRE} = \frac{B8 - B5}{B8 + B5}
$$

- **What it measures:** chlorophyll concentration in already-green
  canopies. The red-edge band B5 sits in the steep slope between
  red-absorption and NIR-reflectance, where the reflectance is very
  sensitive to chlorophyll content.
- **Why we care about it on top of NDVI:** NDRE keeps responding when
  NDVI is saturated, so it is the better metric for telling a *thriving*
  plantation from a merely *green* one. It is also more sensitive to
  early nitrogen stress.
- **Range / interpretation:** typical values 0.1–0.5; higher means more
  chlorophyll per unit canopy.

### 3.3 NDMI — Normalised Difference Moisture Index

$$
\mathrm{NDMI} = \frac{B8 - B11}{B8 + B11}
$$

- **What it measures:** canopy water content. SWIR1 (B11) is strongly
  absorbed by liquid water in leaves; NIR is not. The ratio therefore
  tracks how hydrated the canopy is.
- **Interpretation:** rising NDMI usually means a denser, more turgid
  canopy. Falling NDMI in summer is a drought-stress signal.
- **Why it is in the "Good" rule:** a true establishing forest puts on
  leaf mass *and* water-bearing tissue. A field that just got grass over
  it will green up on NDVI but barely move on NDMI.

### 3.4 BSI — Bare Soil Index

$$
\mathrm{BSI} = \frac{(B11 + B4) - (B8 + B2)}{(B11 + B4) + (B8 + B2)}
$$

- **What it measures:** exposed soil and dry mineral surfaces. Soil is
  bright in SWIR + red and comparatively dim in NIR + blue, so the
  combination above is positive on bare ground and negative under
  vegetation.
- **Range:** typically −0.3 (dense canopy) to +0.3 (bare soil) under
  normal conditions; values around +0.4 to +0.6 are common over fire
  scars or quarries.
- **Why it is in the rules:** a successful establishment must show
  *less* bare ground after a few years. NDVI rising while BSI is also
  rising is a red flag (often grass replacing canopy, or a road / firebreak
  being widened).

### 3.5 NDWI (McFeeters)

$$
\mathrm{NDWI} = \frac{B3 - B8}{B3 + B8}
$$

- **What it measures:** *open water bodies* (lakes, reservoirs, flooded
  fields). Water reflects more green than NIR; vegetation does the
  opposite.
- **Caveat — this is not "vegetation water":** we use the McFeeters
  formulation, which targets standing water. It is intentionally
  *different* from NDMI (which targets leaf water using SWIR). NDWI is
  used in the rules to make sure we do not call a newly flooded paddy
  field "good establishment".
- **What "high NDWI" means on this map:** > ~0.2 → almost certainly open
  water.

### 3.6 NBR — Normalised Burn Ratio

$$
\mathrm{NBR} = \frac{B8 - B12}{B8 + B12}
$$

- **What it measures:** burn severity and post-fire recovery. SWIR2
  (B12) shoots up over burned vegetation (charcoal, dry residue), while
  NIR drops. NBR therefore plummets right after a fire.
- **The change form $\Delta$NBR = NBR_recent − NBR_early** is one of
  the canonical fire-recovery metrics: large positive values mean a
  burned pixel has revegetated; values near zero mean it has not.
- **Why we show it:** on the Portugal AOI most of the dramatic recent
  signal comes from re-greening of the 2017 Pedrógão Grande fire scar.
  NBR_change makes that recovery the most visually obvious layer on
  the map.

---

## 4. The single-period layers

When you open the map, the layers visible by default are the recent
summer composite:

| Layer | What you see | Color ramp |
|---|---|---|
| **Natural color (RGB)** | What a human eye would see from orbit | true colors |
| **NDVI** | Vegetation density now | brown → yellow → green |
| **LAI-e** | Leaf area index estimate (proxy) | white → dark green |
| **FCOVER** | Fractional ground cover (proxy) | white → dark green |
| **BSI** | Bare-soil index now | green → yellow → brown |
| **NDWI** | Open water now | white → blue |
| **NBR** | Burn / recovery signal now | red → yellow → green |

**LAI-e and FCOVER** are linear proxies used in this build (a stand-in
for the full ESA biophysical neural network, which is currently not
wired up). Treat their *spatial pattern* as informative, but do not
read absolute values from them.

The early-period RGB and recent-period RGB layers are hidden by default;
toggle them on to do a quick visual "before/after" without numbers.

---

## 5. The change layers

For each spectral index we compute:

$$
\Delta X = X_{\text{recent}} - X_{\text{early}}
$$

i.e. a per-pixel **subtraction** of the early composite from the recent
composite. Positive $\Delta$ means the index went up; negative means it
went down.

The color ramp is always a **diverging palette** centred on zero, so a
pixel that did not change is rendered white / neutral. The default
range is −0.3 → 0 → +0.3 (clipped at the ends).

| Layer | Palette | Reading |
|---|---|---|
| **NDVI_change** | red → white → green | green = greened up |
| **NDRE_change** | red → white → green | green = chlorophyll up |
| **NDMI_change** | red → white → green | green = canopy got wetter |
| **BSI_change**  | green → white → red | red = more bare soil now (worse) |
| **NDWI_change** | brown → white → blue | blue = more open water; brown = drier |
| **NBR_change**  | red → white → green | green = fire recovery / regrowth |

**Important orientation note for BSI:** because more bare soil is a
*bad* sign, its diverging palette is intentionally inverted relative to
NDVI. Red on the BSI change layer is therefore *always* the bad
direction, just like red on the NDVI change layer is the bad direction.
This consistency was a deliberate design choice — once you internalise
it, you can read every change layer without thinking about the sign.

These layers are **hidden by default**; tick them in the layer panel
to inspect individual signals.

---

## 6. The establishment status map (the headline layer)

This is the layer visible by default in the bottom-right legend. It is
the synthesis of everything above — a single categorical raster with
four classes.

### 6.1 The rule logic

Let $T_{\Delta} = 0.05$, $T_{\text{NDVI low}} = 0.30$,
$T_{\text{BSI high}} = 0.20$, $T_{\text{water}} = 0.20$
(all configurable in the YAML under `classification:`).

```
Good       = (NDVI_change >  T_change)
         AND (BSI_change  < -T_change)
         AND (NDMI_change > -T_change)     # i.e. canopy water did NOT drop
         AND (NDWI_recent <  T_water)      # not standing water

Weak       = (NDVI_recent <  T_ndvi_low)
          OR (BSI_recent  >  T_bsi_high)

Moderate   = any valid pixel that is neither Good nor Weak (default)

Uncertain  = any input band masked over this pixel (no data)
```

Resolution / priority order if a pixel matches more than one rule:

> **Uncertain > Weak > Good > Moderate**

(Implemented in `build_establishment_status()` by starting every pixel
at `Moderate`, promoting matches to `Good`, overriding with `Weak`,
and finally masking to `Uncertain` wherever any input band was missing.)

### 6.2 The palette

| Code | Class | Color (hex) | Meaning |
|---|---|---|---|
| 0 | Uncertain | `#888888` grey | At least one required band was masked (clouds, edge, etc.) |
| 1 | Good establishment | `#1f5e1f` dark green | Greener, drier soil signature, no sign of stress or flooding |
| 2 | Moderate establishment | `#a8d666` light green | Valid pixel, neither thriving nor obviously failing |
| 3 | Weak / bare-soil risk | `#c47a1e` orange | Stayed sparse (NDVI low) or stayed bare (BSI high) — possible failure |

The legend in the bottom-right of the map shows exactly these four
swatches and labels.

### 6.3 How to read it

- **Dark-green patches** are pixels where every component of the rule
  passed: vegetation increased, soil exposure decreased, canopy moisture
  did not collapse, and there is no standing water. That is the closest
  thing the data can give you to "this looks like a successful
  establishment".
- **Light green** is the default fallback. Treat it as "nothing alarming,
  but I cannot confidently say it is thriving" — most stable existing
  forest will end up here because there is no large change between
  early and recent.
- **Orange** is the actionable class. It says either "current NDVI is
  still below 0.3" (i.e. the canopy never closed) or "BSI is above 0.2"
  (i.e. soil is still exposed). On a replanting AOI this is where you
  would want to send someone with boots.
- **Grey** means *we cannot answer for this pixel* — usually persistent
  cloud cover or pixels right on the AOI edge. Do not interpret it as
  bad.

### 6.4 The histogram printed in the console

Every run prints, e.g.

```
Establishment status histogram (portugal-aoi):
  Uncertain                            0 px  (  0.0%)
  Good establishment              55,965 px  ( 51.3%)
  Moderate establishment          49,212 px  ( 45.1%)
  Weak / bare-soil risk            3,858 px  (  3.5%)
```

This is the **areal share of each class** over the full AOI. It is
useful in two ways:

1. As a sanity check that the thresholds are not pathological (e.g. if
   `Weak` ever explodes past 60 %, your thresholds are too strict).
2. As a single-number summary you can quote: "51 % of the AOI shows a
   greening signal consistent with successful establishment".

---

## 7. What this map is and is not

### What it *is*

- A **per-pixel summary of optical change** between two seasonal
  composites built from cloud-masked Sentinel-2 imagery.
- A **transparent rule-based classification** — every class boundary
  is a one-line inequality you can re-read in `build_establishment_status`.
- **Reproducible**: change the YAML, re-run with `--map-only`, get a new
  map in about 30 seconds.

### What it is *not*

- **It is not proof of tree survival.** Optical sensors see canopy
  reflectance. A field of fast-growing grass can mimic the spectral
  signature of a young plantation. Use this map to *prioritise* field
  visits, not to *replace* them.
- **It is not species-specific.** NDVI cannot tell pine from
  eucalyptus from bracken.
- **It is not stress-attribution.** A pixel classified as `Weak` could
  be that way because of drought, pest, fire, browsing, herbicide drift,
  or just slow species. Disambiguating those is the role of the future
  stress-risk layer (Milestone 4 in `vision.md`).
- **It is not calibrated against ground truth.** The thresholds in
  `classification:` are sensible defaults, not measured constants. Once
  you have field reference points they should be tuned per region.

---

## 8. Layer panel — quick reference

When you tick layers on/off in the top-right control:

- **Always look at**: `Establishment status` (the headline product).
- **Cross-check the "Good" patches** by toggling on `NDVI_change`
  (should be green there) and `BSI_change` (should be green there too,
  i.e. less bare soil).
- **Investigate "Weak" patches** by toggling on the recent NDVI layer
  and recent BSI layer — both should be in their "warning" colors there.
- **Sanity-check grey ("Uncertain")** by toggling on the early and
  recent RGB composites. Persistent cloud cover on either period is
  the usual cause.
- **AOI outline** is the red polygon; it is purely cosmetic and has no
  data attached.

---

## 9. What changes if you swap regions

The same map structure is regenerated for any region by changing the
`--config` flag:

| Region | Config | AOI | Reference event |
|---|---|---|---|
| Pedrógão Grande, Portugal | `config_pt.yaml` | ~816 km² | June 2017 fire |
| Rugballegård Skov, Denmark | `config_DK.yaml` | ~6 km² | April 2022 planting |

The Danish AOI is two orders of magnitude smaller, so the same
classification thresholds will look much grainier — but the
interpretation of each color is identical.

---

## 10. Where the math lives in the code

If you want to verify any claim in this document against the
implementation:

| Concept | Function | File |
|---|---|---|
| Cloud Score Plus mask | `build_s2_col` | `src/gee.py` |
| All six spectral indices | `add_spectral_indices` | `src/gee.py` |
| Multi-year seasonal median | `make_seasonal_composite` | `src/gee.py` |
| Establishment rules | `build_establishment_status` | `src/gee.py` |
| Thresholds | `classification:` block | `data/config_*.yaml` |
| Palettes & legend | `ESTAB_CLASSES`, `VIS_CHANGE`, `_legend_html` | `src/gee.py` |
