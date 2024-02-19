def Planet_Downscale_Landsat(img_high, 
                             img_low,
                             resolution_low, resolution_high,
                             window_size, stride):
    '''
    This function is only used for downscaling the Planet-image scale (3m) to Landsat-image scale (30m).
    The image contains four bands.
    
    parameters:
    img_high: Planet image at 3-m resolution with four bands.
    img_low: Downscaled Planet image at 30-m resolution with four bands.
    resolution_low: "Landsat" pixel resolution, 30 m.
    resolution_high: Planet pixel resolution, 3 m.
    window_size: "Landsat" pixel resolution, 30 m.
    stride: The ratio between Landsat and Planet pixel resolution, 30/3=10.
    
    return:
    The downscaled Planet image at 30-m resolution with four bands will be generated.
    '''
    # Necessary libraries
    import numpy as np
    import pandas as pd
    import rasterio
    import os
    import math 
    import time

    # Parameters setting
    scaling_factor = resolution_low/resolution_high
    window_size = int(window_size)
    stride = int(stride)
    
    tic = time.perf_counter()
    with rasterio.open(img_high) as src:
        # Read the image data
        img_og = src.read()
        # Get the geotransform and metadata
        img_og_transform = src.transform
        img_og_metadata = src.meta.copy()
    # Create an empty array for saving the Median value
    upscaled_data_ag1 = np.zeros((img_og.shape[0], img_og.shape[1] // stride, img_og.shape[2] // stride), dtype=img_og.dtype)

    for i in range(0, img_og.shape[1] - window_size + 1, stride):
        for j in range(0, img_og.shape[2] - window_size + 1, stride):
            # Get the local array
            window_r = img_og[0, i:i+window_size, j:j+window_size]
            window_g = img_og[1, i:i+window_size, j:j+window_size]
            window_b = img_og[2, i:i+window_size, j:j+window_size]
            window_nir = img_og[3, i:i+window_size, j:j+window_size]
            # Get the median value
            median_r = np.median(window_r.reshape(-1,1),axis=0)
            median_g = np.median(window_g.reshape(-1,1),axis=0)
            median_b = np.median(window_b.reshape(-1,1),axis=0)
            median_nir = np.median(window_nir.reshape(-1,1),axis=0)
            # Attach to the empty array
            upscaled_data_ag1[0,int(i/10),int(j/10)] = median_r
            upscaled_data_ag1[1,int(i/10),int(j/10)] = median_g
            upscaled_data_ag1[2,int(i/10),int(j/10)] = median_b            
            upscaled_data_ag1[3,int(i/10),int(j/10)] = median_nir

    # Update the geotransform and metadata for the resized image
    resized_transform = (
        resolution_low,
        img_og_transform[1] / scaling_factor,
        img_og_transform[2],
        img_og_transform[3],
        -resolution_low,
        img_og_transform[5])
    resized_metadata = img_og_metadata.copy()
    resized_metadata.update({
        'count': int(img_og.shape[0]),
        'width': int(src.width / scaling_factor),
        'height': int(src.height / scaling_factor),
        'transform': resized_transform})
    # Save the resized image
    with rasterio.open(img_low, 'w', **resized_metadata) as dst:
        dst.write(upscaled_data_ag1)

    toc = time.perf_counter()
    print(f"Processing these files costs {toc - tic:0.4f} seconds")
    print("Downscaling is finished!!!") 
    
    return()