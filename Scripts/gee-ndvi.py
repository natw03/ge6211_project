import ee
import os
import requests

ee.Initialize(project='YOURGEECODEHERE')

project_areas = [
    {"id": 674,  "area": "kalimantan", "geom": ee.Geometry.Rectangle([112.10633568, -3.36245172, 112.50615453, -2.73042218])},
    {"id": 1477, "area": "kalimantan", "geom": ee.Geometry.Rectangle([112.97576164, -3.07360963, 113.34959424, -2.01969795])},
    {"id": 1498, "area": "sumatera",   "geom": ee.Geometry.Rectangle([102.30205101,  0.58183914, 102.68333178,  0.81302914])},
    {"id": 2403, "area": "sumatera",   "geom": ee.Geometry.Rectangle([102.33963344,  0.25842711, 103.05276502,  0.65248461])},
]

YEARS = list(range(2001, 2021)) + [2022, 2023]

output_base = "YOURPATHHERE"

def download_image(image, geom, output_path, scale=1000):
    """Download a GEE image as GeoTIFF to output_path."""
    url = image.getDownloadURL({
        "scale": scale,
        "crs": "EPSG:4326",
        "region": geom,
        "format": "GEO_TIFF",
    })
    response = requests.get(url, stream=True)
    if response.status_code == 200:
        with open(output_path, "wb") as f:
            f.write(response.content)
        print(f"  Saved: {output_path}")
    else:
        print(f"  Failed ({response.status_code}): {output_path}")


def fetch_ndvi(years=YEARS):
    out_dir = os.path.join(output_base, "ndvi")
    os.makedirs(out_dir, exist_ok=True)

    for site in project_areas:
        for year in years:
            print(f"NDVI {site['id']} {year}...")
            image = (ee.ImageCollection("MODIS/061/MOD13A2")
                     .filter(ee.Filter.date(f"{year}-01-01", f"{year}-12-31"))
                     .filterBounds(site["geom"])
                     .select("NDVI")
                     .mean()
                     .multiply(0.0001)
                     .clip(site["geom"]))
            out_path = os.path.join(out_dir, f"NDVI_{site['id']}_{year}.tif")
            download_image(image, site["geom"], out_path, scale=1000)
