def CanopyHeight_Vine(dir_LAI, dir_DEM, dir_R, dir_NIR, 
                      NoDataValue, Veg_threshold, Soil_threshold, height_threshold,
                      dir_output, output_name):
    '''
    LAI image, DEM data, Red band image, Near-infrared band image, and the threshold for NDVI to identify the vegetation pixel. Additionally, images need to be aligned with each other, and the resolution for spectral and DEM need to be the same.<br>
    (1) when there is no canopy pixel, the height is 0;<br>
    (2) when there is canopy pixel:<br>
    &emsp;1) the average height is smaller than 1.4 m, the height is 1.4 m;<br>
    &emsp;2) the average height is bigger than 2.8 m, the height is 2.8 m;<br>
    &emsp;3) the average height is between 1.4 to 2.8, the height is the average value.
    '''
    Array_LAI = arcpy.RasterToNumPyArray(dir_LAI, nodata_to_value=NoDataValue)
    Array_DEM = arcpy.RasterToNumPyArray(dir_DEM, nodata_to_value=NoDataValue)
    Array_R = arcpy.RasterToNumPyArray(dir_R, nodata_to_value=NoDataValue)
    Array_NIR = arcpy.RasterToNumPyArray(dir_NIR, nodata_to_value=NoDataValue)
    Array_NDVI = (Array_NIR-Array_R)/(Array_NIR+Array_R)

    dims_LAI = Array_LAI.shape
    print("LAI dimension is:",dims_LAI[0],dims_LAI[1])
    dims_DEM = Array_DEM.shape
    print("The dimesnion of the high-resolution data is:",dims_DEM[0],dims_DEM[1])
    hor_pixel = int(dims_DEM[0]/dims_LAI[0])
    ver_pixel = int(dims_DEM[1]/dims_LAI[1])
    print("One LAI pixel contains",hor_pixel,"pixels of DEM (also spectral)!")

    # Get the information from LAI map for data output
    fid=gdal.Open(dir_LAI)
    input_lai=fid.GetRasterBand(1).ReadAsArray()
    dims_lai=input_lai.shape
    # Read the GDAL GeoTransform to get the pixel size
    lai_geo=fid.GetGeoTransform()
    lai_prj=fid.GetProjection()
    fid=None
    # Compute the dimensions of the output file
    geo_out=list(lai_geo)
    geo_out=tuple(geo_out)

    # Main part for canopy height calculation
    CanopyHeight = np.empty((dims_LAI[0],dims_LAI[1]))
    CanopyHeight[:] = np.nan
    local_lowest_DEM = 0

    for irow in range(dims_LAI[0]):
        start_row = irow*ver_pixel
        end_row = start_row + hor_pixel
        for icol in range(dims_LAI[1]):
            start_col = icol*hor_pixel
            end_col = start_col + hor_pixel

            local_DEM_new = []
            local_NDVI = Array_NDVI[start_row:end_row,start_col:end_col]
            local_NDVI[local_NDVI < 0] = 0
            local_DEM = Array_DEM[start_row:end_row,start_col:end_col]
            # only soil pixel in the domain, then the canopy height is 0
            # and then renew the lowest height
            if np.nanmax(local_NDVI) < Soil_threshold:
                CanopyHeight[irow,icol] = 0
                local_lowest_DEM = np.nanmin(local_DEM)
            # pixel with NDVI bigger than the threshold for soil NDVI
            else:
                # if the maximum height is smaller than the threshold, make it equal to the threshold
                # if the maximum height is bigger than the threshold, calculate the average value for those higher heights
                # index for soil pixel and non-soil pixel
                index_soil = np.where(local_NDVI <= Soil_threshold)
                index_soil_reverse = np.where(local_NDVI > Soil_threshold)
                index_veg = np.where(local_NDVI > Veg_threshold)
                # if there are soil pixel, renew the height and then provide corresponding parameters
                if len(index_soil[0]) > 0:
                    # renew the height
                    local_lowest_DEM = np.nanmin(local_DEM)
                    local_DEM_new = local_DEM - local_lowest_DEM
                    local_DEM_new[local_DEM_new < height_threshold] = np.nan
                    if np.nanmax(local_NDVI) < Veg_threshold:
                        CanopyHeight[irow,icol] = 0                    
                    else:
                        if len(index_veg[0])/(hor_pixel*ver_pixel) <= 0.05:
                            CanopyHeight[irow,icol] = 0
                        else:
                            if np.isnan(local_DEM_new).all():
                                CanopyHeight[irow,icol] = height_threshold
                            else:
                                CanopyHeight[irow,icol] = np.nanmean(local_DEM_new)
                # if there are no soil pixel, using the previous height and then provide corresponding parameters
                else:
                    local_DEM_new = local_DEM - local_lowest_DEM
                    local_DEM_new[local_DEM_new < height_threshold] = np.nan
                    if np.nanmax(local_NDVI) < Veg_threshold:
                        CanopyHeight[irow,icol] = 0                    
                    else:
                        if len(index_veg[0])/(hor_pixel*ver_pixel) <= 0.05:
                            CanopyHeight[irow,icol] = 0
                        else:
                            if np.isnan(local_DEM_new).all():
                                CanopyHeight[irow,icol] = height_threshold
                            else:
                                CanopyHeight[irow,icol] = np.nanmean(local_DEM_new)

    plt.figure()
    plt.scatter(tmp_hc,tmp_hc)
    plt.title("Check the canopy-height range from the 1:1 line.")
    plt.show()

    # Write the output file
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(dir_output+"\\"+output_name, dims_LAI[1], dims_LAI[0], 1, gdal.GDT_Float32)
    ds.SetGeoTransform(geo_out)
    ds.SetProjection(lai_prj)
    band=ds.GetRasterBand(1)
    band.WriteArray(CanopyHeight)
    band.SetNoDataValue(NoDataValue)
    band.FlushCache()
    ds = None
    return("Canopy-height calculation is finished.")