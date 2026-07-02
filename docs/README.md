# Interactive network viewer

`index.html` is a self-contained interactive viewer for the hazard-chain
networks (uses [`vis-network`](https://visjs.github.io/vis-network/) from a CDN).

- **Figure 1** — marine chain, fine taxonomy (39 nodes / 152 edges)
- **Figure 2** — marine (15 classes / 57 edges) vs mainland China (15 / 47)
- Click a node to light up its neighbours and dim the rest.
- Colour by four-layer taxonomy or by degree; Chinese / English toggle.

`graph-data.js` is generated from the `data/` cascade tables; the edge sets
reproduce the manuscript exactly.

Served via GitHub Pages (Settings -> Pages -> Deploy from `main` / `docs`).
