
# Re-check sensitivity with new thresholds
```
python src/calibrate_thresholds.py --config config_pt.yaml --factors "0.9,1.0,1.1" --csv
```

output: 
- updates sensitivity table and re-writes the sensitivity_portugal-aoi.csv

# Rebuild map with new thresholds

```
python src/gee.py --config config_pt.yaml --map-only
```

output : the html map and prints the new class histograms

# Full run with time series plots

```
python src/gee.py --config config_pt.yaml
```

