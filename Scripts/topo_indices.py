import os
import glob
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.merge import merge
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling

base_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(base_dir, '..', 'data')
sumatera_tiles_dir = os.path.join(data_dir, 'indices', 'DEM', 'sumatera_tiles')
kalimantan_tiles_dir = os.path.join(data_dir, 'indices', 'DEM', 'kalimantan_tiles') 
study_areas_dir = os.path.join(data_dir, 'indices_boundary')

regions = {
    'sumatera': {
        'tiles_dir': sumatera_tiles_dir,
        'shp_path': os.path.join(study_areas_dir, 'bound_peat_sumatera.shp'),
        'epsg': 'EPSG:23831',
        'out_dem': os.path.join(data_dir, 'indices', 'DEM', 'dem_sumatera.tif'),
        'out_slope': os.path.join(data_dir, 'indices', 'Slope', 'slope_sumatera.tif'),
        'out_aspect': os.path.join(data_dir, 'indices', 'Aspect', 'aspect_sumatera.tif')
    },
    'kalimantan': {
        'tiles_dir': kalimantan_tiles_dir,
        'shp_path': os.path.join(study_areas_dir, 'bound_peat_kali.shp'),
        'epsg': 'EPSG:23832',
        'out_dem': os.path.join(data_dir, 'indices', 'DEM','dem_kali.tif'),
        'out_slope': os.path.join(data_dir,'indices', 'Slope','slope_kali.tif'),
        'out_aspect': os.path.join(data_dir, 'indices', 'Aspect','aspect_kali.tif')
    }
}

def process_terrain(region_name, config):
    print(f"\n{'='*40}")
    print(f"Processing Terrain for {region_name.upper()}")
    print(f"{'='*40}")

    search_criteria = os.path.join(config['tiles_dir'], "*.jp2")
    dem_fps = glob.glob(search_criteria)

    if not dem_fps:
        print(f"[!] No DEM tiles found for {region_name} in {config['tiles_dir']}")
        return
    
    print("Loading boundary and aligning coordinate systems...")
    boundary_shp = gpd.read_file(config['shp_path'])

    with rasterio.open(dem_fps[0]) as first_dem:
        dem_crs = first_dem.crs
    boundary_dem_crs = boundary_shp.to_crs(dem_crs)
    bounds = tuple(boundary_dem_crs.total_bounds)

    print(f"Found {len(dem_fps)} DEM tile(s).")
    
    if len(dem_fps) == 1:
        print("Only 1 tile found! Skipping merge overhead and clipping directly...")
        with rasterio.open(dem_fps[0]) as src:
            out_meta = src.meta.copy()
            clipped_image, clipped_transform = mask(
                src, boundary_dem_crs.geometry.values, crop=True, filled=True, nodata=-9999
            )
    else:
        print("Merging tiles strictly within the regional bounding box...")
        src_files_to_mosaic = []
        try:
            for fp in dem_fps:
                src = rasterio.open(fp)
                src_files_to_mosaic.append(src)
                
            mosaic, out_trans = merge(src_files_to_mosaic, bounds=bounds, nodata=-9999)
            
            out_meta = src_files_to_mosaic[0].meta.copy()
            with rasterio.io.MemoryFile() as memfile:
                out_meta.update({"driver": "GTiff", "height": mosaic.shape[1], "width": mosaic.shape[2], "transform": out_trans, "nodata": -9999})
                with memfile.open(**out_meta) as mem_src:
                    mem_src.write(mosaic)
                    clipped_image, clipped_transform = mask(
                        mem_src, boundary_dem_crs.geometry.values, crop=True, filled=True, nodata=-9999
                    )
        finally:
            for src in src_files_to_mosaic:
                src.close()

    out_meta.update({
        "driver": "GTiff", 
        "height": clipped_image.shape[1], 
        "width": clipped_image.shape[2], 
        "transform": clipped_transform, 
        "nodata": -9999
    })

    print(f"Reprojecting clipped DEM to {config['epsg']} (Meters) for accurate math...")
    target_crs = config['epsg']
    transform, width, height = calculate_default_transform(
        dem_crs, target_crs, clipped_image.shape[2], clipped_image.shape[1], *boundary_dem_crs.total_bounds
    )
    
    proj_dem = np.zeros((height, width), dtype='float32')
    reproject(
        source=clipped_image[0], 
        destination=proj_dem, 
        src_transform=clipped_transform,
        src_crs=dem_crs, 
        dst_transform=transform, 
        dst_crs=target_crs,
        resampling=Resampling.bilinear, 
        src_nodata=-9999, 
        dst_nodata=-9999
    )

    final_meta = out_meta.copy()
    final_meta.update({
        "transform": transform, 
        "crs": target_crs, 
        "width": width, 
        "height": height, 
        "dtype": 'float32'
    })
    
    os.makedirs(os.path.dirname(config['out_dem']), exist_ok=True)
    with rasterio.open(config['out_dem'], "w", **final_meta) as dest:
        dest.write(proj_dem, 1)
    print(f"Saved strictly isolated regional DEM to {config['out_dem']}")

    print("Calculating Slope and Aspect...")
    proj_dem[proj_dem == -9999] = np.nan

    dx = transform[0]
    dy = -transform[4]

    y_grad, x_grad = np.gradient(proj_dem, dy, dx)

    slope = np.degrees(np.arctan(np.sqrt(x_grad**2 + y_grad**2)))
    slope[np.isnan(proj_dem)] = -9999

    os.makedirs(os.path.dirname(config['out_slope']), exist_ok=True)
    with rasterio.open(config['out_slope'], "w", **final_meta) as dest:
        dest.write(slope.astype('float32'), 1)
    print(f"Saved Slope to {config['out_slope']}")

    aspect = np.degrees(np.arctan2(-x_grad, y_grad))
    aspect = (aspect + 360) % 360
    aspect[slope < 0.01] = -1 
    aspect[np.isnan(proj_dem)] = -9999

    os.makedirs(os.path.dirname(config['out_aspect']), exist_ok=True)
    with rasterio.open(config['out_aspect'], "w", **final_meta) as dest:
        dest.write(aspect.astype('float32'), 1)
    print(f"Saved Aspect to {config['out_aspect']}")

if __name__ == "__main__":
    for region, conf in regions.items():
        process_terrain(region, conf)
    print("\nAll terrain processing complete")