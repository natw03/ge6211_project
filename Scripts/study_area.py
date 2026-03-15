import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.features import shapes

# --- Setup Paths ---
base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, '..', 'data')
study_areas_dir = os.path.join(data_dir, 'study_areas_updated')

original_raster_path = os.path.join(data_dir, 'peatland.tif')

# Region-specific configurations
regions = {
    'borneo': {
        'shp_path': os.path.join(study_areas_dir, 'study_areas_borneo.shp'),
        'epsg': 'EPSG:23835',
        'out_peat': os.path.join(data_dir, 'bound_peat_borneo.shp'),
        'out_buffer': os.path.join(data_dir, 'buffer_5km_borneo.shp')
    },
    'sumatra': {
        'shp_path': os.path.join(study_areas_dir, 'study_areas_sumatra.shp'),
        'epsg': 'EPSG:23831',
        'out_peat': os.path.join(data_dir, 'bound_peat_sumatra.shp'),
        'out_buffer': os.path.join(data_dir, 'buffer_5km_sumatra.shp')
    }
}

def process_region(region_name, config):
    print(f"\n{'='*40}")
    print(f"--- Processing {region_name.upper()} ---")
    print(f"{'='*40}")
    
    target_crs = config['epsg']
    
    # 1. Load Region-Specific Shapefile
    print(f"Loading {region_name} shapefile and projecting to {target_crs}...")
    region_shp = gpd.read_file(config['shp_path'])
    region_proj = region_shp.to_crs(target_crs)
    
    # 2. Create the 5km Buffer
    print("Creating 5km buffer...")
    region_buffer = region_proj.copy()
    region_buffer['geometry'] = region_proj.geometry.buffer(5000)
    
    try:
        region_buffer.to_file(config['out_buffer'], mode='w')
    except PermissionError:
        print(f"[!] Warning: Could not save buffer. Close {config['out_buffer']} if open.")

    # 3. Apply the 5km Buffer as a PIXEL MASK
    print("Masking master raster using the exact 5km buffer geometries...")
    buffer_wgs84 = region_buffer.to_crs('EPSG:4326')
    
    with rasterio.open(original_raster_path) as src:
        out_image, out_transform = mask(
            src, 
            buffer_wgs84.geometry.values, 
            crop=True, 
            nodata=0 
        )
        
        profile = src.profile.copy()
        profile.update({
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
            "nodata": 0
        })
        
        with rasterio.io.MemoryFile() as memfile:
            with memfile.open(**profile) as mem_src:
                mem_src.write(out_image)
                
                # 4. Reproject raster to local CRS
                print(f"Reprojecting masked raster to {target_crs}...")
                transform, width, height = calculate_default_transform(
                    mem_src.crs, target_crs, mem_src.width, mem_src.height, *mem_src.bounds
                )
                
                destination = np.zeros((height, width), mem_src.dtypes[0])
                reproject(
                    source=rasterio.band(mem_src, 1),
                    destination=destination,
                    src_transform=mem_src.transform,
                    src_crs=mem_src.crs,
                    dst_transform=transform,
                    dst_crs=target_crs,
                    resampling=Resampling.nearest
                )
                
                # 5. Polygonize
                print("Polygonizing strictly buffered peatland pixels...")
                peat_mask = (destination == 1)
                geoms = list(
                    {'properties': {'peatland': v}, 'geometry': s}
                    for i, (s, v) 
                    in enumerate(shapes(destination, mask=peat_mask, transform=transform)) 
                    if v == 1
                )
                
    if not geoms:
        print(f"No peatland found within the {region_name} 5km buffer.")
        return

    # 6. SLEDGEHAMMER: Strict Vector Clipping and Smoothing
    print("Executing final vector clip to force absolute regional isolation...")
    peat_gdf = gpd.GeoDataFrame.from_features(geoms, crs=target_crs)
    
    # Fix any minor topology errors before clipping
    peat_gdf['geometry'] = peat_gdf.geometry.buffer(0)
    
    # THIS IS THE FIX: Physically chop off anything outside the 5km local buffer
    peat_gdf = gpd.clip(peat_gdf, region_buffer)
    
    print("Dissolving and smoothing polygons...")
    peat_gdf = peat_gdf.dissolve(by='peatland').reset_index()
    peat_gdf['geometry'] = peat_gdf.geometry.simplify(tolerance=100, preserve_topology=True)
    peat_gdf = peat_gdf[~peat_gdf.geometry.is_empty]
    
    # 7. Save
    try:
        peat_gdf.to_file(config['out_peat'], mode='w')
        print(f"Successfully saved smooth, absolutely isolated shapefile to {config['out_peat']}")
    except PermissionError:
        print(f"\n[!] PERMISSION ERROR: Could not overwrite '{config['out_peat']}'.")
        print("[!] Please ensure the file is NOT open in QGIS/ArcGIS.")
        
    print(f"Finished {region_name.capitalize()}.\n")

# --- Execute Script ---
if __name__ == "__main__":
    for region, conf in regions.items():
        process_region(region, conf)
        
    print("All processing complete!")