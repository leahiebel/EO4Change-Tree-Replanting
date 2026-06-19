const map = L.map("map", { preferCanvas: true }).setView([39.7, -8.0], 7);

const statusEl = document.getElementById("status");
const legendEl = document.getElementById("legend");
const opacityEl = document.getElementById("opacity");

const osm = L.tileLayer("https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png", {
  maxZoom: 19,
  attribution: "&copy; OpenStreetMap contributors",
}).addTo(map);

L.control.scale({ metric: true, imperial: false }).addTo(map);

const overlayLayers = {};
const rasterLayers = [];
const layerControl = L.control.layers({ OpenStreetMap: osm }, overlayLayers, { collapsed: false }).addTo(map);

function setStatus(message) {
  statusEl.textContent = message;
}

function renderLegend(layer) {
  legendEl.innerHTML = "";
  const legend = layer.legend || {};
  if (legend.classes) {
    Object.entries(legend.classes).forEach(([value, item]) => {
      if (value === "0") return;
      const row = document.createElement("div");
      row.className = "legend-row";
      row.innerHTML = `<span class="swatch" style="background:${item.color}"></span><span>${item.label}</span>`;
      legendEl.appendChild(row);
    });
    return;
  }

  const palette = legend.palette || [];
  palette.forEach((color, index) => {
    const row = document.createElement("div");
    row.className = "legend-row";
    const label =
      index === 0
        ? `${legend.min} ${layer.unit}`
        : index === palette.length - 1
          ? `${legend.max} ${layer.unit}`
          : "";
    row.innerHTML = `<span class="swatch" style="background:${color}"></span><span>${label}</span>`;
    legendEl.appendChild(row);
  });
}

function addRasterLayer(layer, addToMap = false) {
  const tile = L.tileLayer(layer.tile_url, {
    minZoom: layer.min_zoom || 5,
    maxZoom: layer.max_zoom || 14,
    opacity: layer.opacity || 0.75,
    attribution: layer.attribution || "",
  });
  tile.on("loading", () => setStatus(`Loading ${layer.title}...`));
  tile.on("load", () => setStatus("Layers ready."));
  tile.on("tileerror", () => setStatus(`Error loading ${layer.title}. Check tile service and COG URL.`));
  tile._eo4c = layer;
  rasterLayers.push(tile);
  layerControl.addOverlay(tile, layer.title);
  if (addToMap) {
    tile.addTo(map);
    renderLegend(layer);
  }
}

function addBoundary(url, name, style) {
  fetch(url)
    .then((response) => {
      if (!response.ok) throw new Error(`${url} not available`);
      return response.json();
    })
    .then((geojson) => {
      const boundary = L.geoJSON(geojson, {
        style,
        onEachFeature: (feature, layer) => {
          const props = feature.properties || {};
          const label = props.site_name || props.project_id || props.feature_id || name;
          layer.bindPopup(`<strong>${name}</strong><br>${label}`);
        },
      }).addTo(map);
      layerControl.addOverlay(boundary, name);
      if (name === "AOI boundary") map.fitBounds(boundary.getBounds(), { padding: [20, 20] });
    })
    .catch(() => setStatus(`${name} boundary not loaded.`));
}

function queryPoint(layer, latlng) {
  if (!layer.cog_url || !layer.tile_url.includes("tile-server.example.com") === false) {
    return Promise.resolve(null);
  }
  const base = layer.tile_url.split("/tiles/")[0];
  const url = `${base}/point/${latlng.lng},${latlng.lat}?url=${encodeURIComponent(layer.cog_url)}`;
  return fetch(url)
    .then((response) => (response.ok ? response.json() : null))
    .catch(() => null);
}

map.on("overlayadd", (event) => {
  if (event.layer && event.layer._eo4c) renderLegend(event.layer._eo4c);
});

map.on("click", (event) => {
  const active = rasterLayers.find((layer) => map.hasLayer(layer));
  if (!active) return;
  const meta = active._eo4c;
  queryPoint(meta, event.latlng).then((payload) => {
    const value = payload && payload.values ? payload.values[0] : "Value query not supported by this tile endpoint";
    L.popup()
      .setLatLng(event.latlng)
      .setContent(`<strong>${meta.title}</strong><br>${value} ${meta.unit || ""}`)
      .openOn(map);
  });
});

opacityEl.addEventListener("input", (event) => {
  const opacity = Number(event.target.value);
  rasterLayers.forEach((layer) => layer.setOpacity(opacity));
});

Promise.all([
  fetch("../data/outputs/layer_manifest.json").then((response) => {
    if (!response.ok) throw new Error("layer_manifest.json not found");
    return response.json();
  }),
])
  .then(([manifest]) => {
    manifest.forEach((layer, index) => addRasterLayer(layer, index === 0));
    addBoundary("aoi.geojson", "AOI boundary", { color: "#111111", weight: 2, fillOpacity: 0 });
    addBoundary("forest.geojson", "Forest polygons", { color: "#238b45", weight: 2, fillOpacity: 0.08 });
    setStatus("Layers ready.");
  })
  .catch((error) => {
    setStatus(`Could not load layer manifest: ${error.message}`);
  });

