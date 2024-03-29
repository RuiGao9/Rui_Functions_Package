def TemperatureSeparation(dir_LAI, dir_RGBNIR, dir_Tr, dir_DSM, 
                          NoDataValue, Veg_threshold, Soil_threshold,
                          band_R, band_NIR, 
                          cellsize_resample, 
                          Azimuth, Altitude, 
                          MiddleProducts="No"):
    '''
    This function is trying to consider the effect of shadow pixel. But the algorithm used in this function cannot identify the shadow pixel since
    the upscaled DSM data (from 0.15 meter pixel to 0.6 meter pixel) cannot tell any shadow pixel. In other words, shadow information is smoothed 
    due to the DSM data upscaling.
    
    Parameters used in this function:
    dir_LAI: the file path of the LAI image, resolution is 3.6 meter by 3.6 meter.
    dir_RGBNIR: the file path of the optical image containing R, G, B, and NIR bands, and the resolution is 0.15 meter b 0.15 meter.
    dir_Tr: the file path of the temperature image in unith of degree C, and the resolution is 0.6 meter by 0.6 meter.
    dir_DSM: the file path of the DSM image in meter and the resolution is 0.15 meter by 0.15 meter.
    NoDataValue: assigne a value to represent the NaN value.
    Veg_threshold: any NDVI pixel value above this threshold represents vegetation pixels.
    Soil_threshold: any NDVI pixel value below this threshold represents soil pixels.
    band_R: the number of layer in the optical image (multiple bands) representing the Red band.
    band_NIR: the number of layer in the optical image (multiple bands) representing the Near-infrared band.
    cellsize_resample: 0.6 meter by 0.6 meter resolution in order to calculate the vine shadow and temperature separation.
    Azimuth: a parameter used for vine shadow calculation.
    Altitude: a parameter used for vine shadow calculation.
    MiddleProducts: default is "No", which means to delete the middle products. Other parameters, like "Yes" will save the middle products.
    '''
    # import libraries
    import arcpy
    import gdal
    import os
    import numpy as np
    import pandas as pd
    from scipy.stats import linregress
    import matplotlib.pyplot as plt
    
    # Optical image processing
    # resample the optical image
    [res_x,res_y] = TellResolution(dir_RGBNIR)
    name_resample = "resample_RGBNIR.tif"
    arcpy.Resample_management(in_raster=dir_RGBNIR, 
                              out_raster=dir_output+"\\"+name_resample,
                              cell_size=str(cellsize_resample)+" "+str(cellsize_resample), 
                              resampling_type="BILINEAR")
    # clip the optical and thermal image the same as the LAI
    extent = TellExtent(dir_LAI)
    name_clip_opt = "clip_RGBNIR.tif"
    arcpy.Clip_management(in_raster=dir_output+"\\"+name_resample, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+name_clip_opt, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")

    # Temperature data clip
    # Resolution of the temperature is 0.6 meter by 0.6 meter
    # The unit is in degree C
    name_clip_tr = "clip_tr.tif"
    arcpy.Clip_management(in_raster=dir_Tr, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+name_clip_tr, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")

    # DSM data clip
    # resample the optical image
    [res_x,res_y] = TellResolution(dir_DSM)
    name_resample_dsm = "resample_DSM.tif"
    arcpy.Resample_management(in_raster=dir_DSM, 
                              out_raster=dir_output+"\\"+name_resample_dsm,
                              cell_size=str(cellsize_resample)+" "+str(cellsize_resample), 
                              resampling_type="BILINEAR")
    # Clip the DSM data
    name_clip_dsm = "clip_DSM.tif"
    arcpy.Clip_management(in_raster=dir_output+"\\"+name_resample_dsm, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+name_clip_dsm, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")

    # Temperature separation
    Array_LAI = arcpy.RasterToNumPyArray(dir_LAI, nodata_to_value=NoDataValue)
    Array_RGBNIR = arcpy.RasterToNumPyArray(dir_output+"\\"+name_clip_opt, nodata_to_value=NoDataValue)
    Array_R = Array_RGBNIR[band_R,:,:]
    Array_NIR = Array_RGBNIR[band_NIR,:,:]

    Array_NDVI = (Array_NIR-Array_R)/(Array_NIR+Array_R)
    Array_NDVI[Array_NDVI<0] = np.nan
    Array_NDVI[Array_NDVI>1] = np.nan
    Array_Tr = arcpy.RasterToNumPyArray(dir_output+"\\"+name_clip_tr, nodata_to_value=NoDataValue)
    Array_Tr[Array_Tr<0] = np.nan
#     print(Array_LAI.shape, Array_RGBNIR.shape, Array_R.shape, Array_NIR.shape, Array_Tr.shape)
#     print(np.nanmax(Array_NDVI),np.nanmin(Array_NDVI))

    # Stefan-Boltzmann Law
    Array_Tr = Array_Tr ** 4

    dims_LAI = Array_LAI.shape
    print("Column of the LAI:",dims_LAI[0],"Row of the LAI:",dims_LAI[1])
    dims_NDVI = Array_NDVI.shape
    print("Column of the spectral data:",dims_NDVI[0],"Row of the spectral data:",dims_NDVI[1])
    hor_pixel = int(dims_NDVI[0]/dims_LAI[0])
    ver_pixel = int(dims_NDVI[1]/dims_LAI[1])
    print("Each LAI pixel contains",hor_pixel,"(column/column) by",ver_pixel,"(row/row) pixels.")

    # Hillshade calculation
    # Hillshade calculation based on the cliped RGBNIR data
    name_hillshade = "Hillshade.tif"
    arcpy.gp.HillShade_sa(dir_output+"\\"+name_clip_dsm, 
                          dir_output+"\\"+name_hillshade, 
                          str(Azimuth), 
                          str(Altitude), 
                          "SHADOWS", "1")
    # Read the hillshade as array to ignore the shadow pixels
    Array_Hillshade = arcpy.RasterToNumPyArray(dir_output+"\\"+name_hillshade, nodata_to_value=NoDataValue)

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

    t_canopy = np.empty((dims_LAI[0],dims_LAI[1]))
    t_canopy[:] = np.nan
    t_soil = np.empty((dims_LAI[0],dims_LAI[1]))
    t_soil[:] = np.nan
    t_coeff = np.empty((dims_LAI[0],dims_LAI[1]))
    t_coeff[:] = np.nan
    print("Dimension of the canopy temperature is:",t_canopy.shape[0],t_canopy.shape[1])
    print("Dimension of the soil temperature is:",t_soil.shape[0],t_soil.shape[1])

    # initial values for these four variables
    renew_slope = NoDataValue
    renew_intercept = NoDataValue
    renew_coeff = NoDataValue
    slope = NoDataValue
    intercept = NoDataValue
    correlation = NoDataValue
    pvalue = NoDataValue
    stderr = NoDataValue

    for irow in range(dims_LAI[0]):
        start_row = irow * hor_pixel
        end_row = start_row + (hor_pixel)
        for icol in range(dims_LAI[1]):
            start_col = icol * ver_pixel
            end_col = start_col + (ver_pixel)

            # Using the hillshade to eliminate the shadow pixel
            local_NDVI = Array_NDVI[start_row:end_row,start_col:end_col]
            local_NDVI[local_NDVI < 0] = NoDataValue
            local_Tr = Array_Tr[start_row:end_row,start_col:end_col]
            local_Hillshade = Array_Hillshade[start_row:end_row,start_col:end_col]
            local_Hillshade[local_Hillshade<=0] = 0
            local_Hillshade[local_Hillshade>0] = 1

            local_NDVI = local_NDVI*local_Hillshade

            tmp_NDVI = local_NDVI.reshape(-1)
            tmp_NDVI[tmp_NDVI<=0] = np.nan
    #         print("Invalid pixel number is:",np.count_nonzero(np.isnan(tmp_NDVI)))
            tmp_Tr = local_Tr.reshape(-1)
            tmp_Tr = np.sqrt(np.sqrt(tmp_Tr)) # unit in degree C

            df = pd.DataFrame()
            df = pd.DataFrame({'NDVI': tmp_NDVI,'Tr': tmp_Tr})
            df = df.dropna()
            df = df.apply(pd.to_numeric, errors='coerce')

            # do regression if valid data existed in the data frame
            if len(df) != 0:
                slope,intercept,correlation,pvalue,stderr = linregress(df['NDVI'],df['Tr'])        
                # slope: slope of the regression
                # intercept: intercept of the regression line
                # correlation: correlation coefficient
                # pvalue: two-sided p-value for a hypothesis test whose null hypothesis is that the slope is zero
                # stderr: standard error of the estimate
            else: pass

            # renew the slope and the intercept if the slope is negative
            if np.nanmean(slope) < 0:
                renew_slope = slope
                renew_intercept = intercept
                renew_coeff = correlation
            else: pass

            # gain index for soil and canopy pixel for each local domain
            index_soil = np.where(local_NDVI <= Soil_threshold)
            index_veg = np.where(local_NDVI >= Veg_threshold)
            # when the domain contains both vegetation and soil
            if len(index_soil[0]) > 0 and len(index_veg[0]) > 0:
                t_canopy[irow,icol] = np.nanmean(local_Tr[index_veg[0],index_veg[1]])
                t_soil[irow,icol] = np.nanmean(local_Tr[index_soil[0],index_soil[1]])
            # when the domain contains vegetation but no soil: estimate the soil temperature
            elif len(index_soil[0]) == 0 and len(index_veg[0]) > 0:
                t_canopy[irow,icol] = np.nanmean(local_Tr[index_veg[0],index_veg[1]])
                t_soil[irow,icol] = ((renew_intercept + renew_slope * Soil_threshold)**2)**2
            # when the domain contains soil but no vegetation: vegetation temperature is "NAN"
            elif len(index_soil[0]) > 0 and len(index_veg[0]) <= 0:
                t_canopy[irow,icol] = np.nan
                t_soil[irow,icol] = np.nanmean(local_Tr[index_soil[0],index_soil[1]])
            # when the domain contains either pure soil or vegetation
            # estimate the soil and vegetation temperature
            elif len(index_soil[0]) == 0 and len(index_veg[0]) == 0:
                t_canopy[irow,icol] = np.nan
                t_soil[irow,icol] = np.nan

            t_coeff[irow,icol] = renew_coeff

    tt_canopy = np.sqrt(np.sqrt(t_canopy.copy())) + 273.15
    tt_soil = np.sqrt(np.sqrt(t_soil.copy())) + 273.15

    # Write the output file
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(dir_output+"\\"+output_name, dims_LAI[1], dims_LAI[0], 3, gdal.GDT_Float32)
    ds.SetGeoTransform(geo_out)
    ds.SetProjection(lai_prj)
    band=ds.GetRasterBand(1)
    band.WriteArray(tt_canopy)
    band.SetNoDataValue(NoDataValue)
    band.FlushCache()
    band=ds.GetRasterBand(2)
    band.WriteArray(tt_soil)
    band.SetNoDataValue(NoDataValue)
    band.FlushCache()
    band=ds.GetRasterBand(3)
    band.WriteArray(t_coeff)
    band.SetNoDataValue(NoDataValue)
    band.FlushCache()
    ds = None
    print("Done!!! Temperature separation is finished.")
    
    # delete the middle products
    if MiddleProducts == "No":
        os.remove(dir_output+"\\"+name_resample)
        os.remove(dir_output+"\\"+name_clip_opt)
        os.remove(dir_output+"\\"+name_clip_tr)
        os.remove(dir_output+"\\"+name_resample_dsm)
        os.remove(dir_output+"\\"+name_clip_dsm)
        os.remove(dir_output+"\\"+name_hillshade)
    else:
        pass
    
    return("Temperature separation is finished!!!")
