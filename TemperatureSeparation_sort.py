def TemperatureSeparation_sort(dir_LAI, dir_RGBNIR, dir_Tr, dir_DSM, 
                               NoDataValue, Veg_threshold, Soil_threshold,
                               band_R, band_NIR, 
                               cellsize_resample, 
                               dir_output, output_name_multiple,
                               Azimuth, Altitude, 
                               MiddleProducts="No"):    
    # import libraries
    import arcpy
    import gdal
    import os
    import numpy as np
    import pandas as pd
    from sklearn import linear_model
    import matplotlib.pyplot as plt


    # RGB-NIR image processing
    # resample the optical image - from 0.15 m to 0.60 m pixel
    [res_x,res_y] = TellResolution(dir_RGBNIR)
    name_resample = "resample_RGBNIR.tif"
    tmp = arcpy.ia.Resample(dir_RGBNIR, "Average", res_x, cellsize_resample)
    tmp.save(dir_output+"\\"+name_resample)
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
    # Temperature image processing
    # The resolution is 0.6 m pixel
    # The unit now is degree C
    name_clip_tr = "clip_tr.tif"
    arcpy.Clip_management(in_raster=dir_Tr, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+name_clip_tr, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")
    # Shadow pixels identifying based on DSM image
    name_hillshade = "Hillshade.tif"
    arcpy.gp.HillShade_sa(dir_DSM, 
                          dir_output+"\\"+name_hillshade, 
                          str(Azimuth), 
                          str(Altitude), 
                          "SHADOWS", "1")
    # Aggregate the Hillshade image,
    # If the 0.6 meter grid contains at least one 0.15 meter shadow pixel,
    # this 0.6 meter grid is shadow pixel
    name_aggregate = "Aggregate.tif"
    arcpy.gp.Aggregate_sa(dir_output+"\\"+name_hillshade, 
                          dir_output+"\\"+name_aggregate, 
                          "4", 
                          "MINIMUM", 
                          "EXPAND", "DATA")
    # 4 is a constant value here, since the resolution for the original DSM and Temperature image is 0.15 and 0.6 meter resolution.
    # Clip the aggregated DSM data
    # This image can be used to identify the shadow pixel: 0 represents the shadow
    name_clip_dsm = "clip_aggregate.tif"
    arcpy.Clip_management(in_raster=dir_output+"\\"+name_aggregate, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+name_clip_dsm, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")
    # Read this image as an array
    # Classify this image containing shadow pixel and non-shadow pixel
    Array_Shadow = arcpy.RasterToNumPyArray(dir_output+"\\"+name_clip_dsm, nodata_to_value=NoDataValue)
    Array_Shadow[Array_Shadow>0] = 1 # non-shadow pixel
    # print(np.nanmax(Array_Shadow),np.nanmin(Array_Shadow))
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
    print(Array_LAI.shape, Array_RGBNIR.shape, Array_R.shape, Array_NIR.shape, Array_Tr.shape)
    print(round(np.nanmax(Array_NDVI),1),round(np.nanmin(Array_NDVI),1))

    # Stefan-Boltzmann Law
    Array_Tr = Array_Tr ** 4

    dims_LAI = Array_LAI.shape
    print("Column of the LAI:",dims_LAI[0],"Row of the LAI:",dims_LAI[1])
    dims_NDVI = Array_NDVI.shape
    print("Column of the spectral data:",dims_NDVI[0],"Row of the spectral data:",dims_NDVI[1])
    hor_pixel = int(dims_NDVI[0]/dims_LAI[0])
    ver_pixel = int(dims_NDVI[1]/dims_LAI[1])
    print("Each LAI pixel contains",hor_pixel,"(column/column) by",ver_pixel,"(row/row) pixels.")

    ## Main part of the temperature separation
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

    [renew_slope,renew_intercept] = [NoDataValue,NoDataValue]
    # initial values for these four variables
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
            local_Shadow = Array_Shadow[start_row:end_row,start_col:end_col]
            num_zero = np.count_nonzero(local_Shadow==0)

            local_NDVI = local_NDVI*local_Shadow

            tmp_NDVI = local_NDVI.reshape(-1)
            tmp_NDVI[tmp_NDVI<=0] = np.nan
            tmp_Tr = local_Tr.reshape(-1)
            tmp_Tr = np.sqrt(np.sqrt(tmp_Tr)) # unit in degree C

            df = pd.DataFrame()
            df = pd.DataFrame({'NDVI': tmp_NDVI,'Tr': tmp_Tr})
            df = df.dropna()
            df = df.apply(pd.to_numeric, errors='coerce')
    #         print(df['NDVI'].values.shape,type(df['NDVI']),type(df['NDVI'].values),type(df[['NDVI']]))
            # Sorting NDVI and Tr
            df = pd.DataFrame({'NDVI': df['NDVI'].sort_values(ascending=True).values,'Tr': df['Tr'].sort_values(ascending=False).values})

            # do regression if valid data existed in the data frame
            if len(df) != 0:
                # Robust linear relationship
                robust_model = linear_model.RANSACRegressor()
                robust_model.fit(df[['NDVI']],df[['Tr']])
                robust_y = robust_model.predict(df[['NDVI']])
                # processing for robust_y
                [num, robust_y_convert] = [len(robust_y),[]]
                for irou in range(0,num):
                    robust_y_convert.append(robust_y[irou][0])
                # correlation coefficient between the temperature predictions and observations
                coef = gfit(df['Tr'].values,np.asarray(robust_y_convert),type_statistic='4',residual='No')
                # getting the slope and intercept from the robust linear regression
                slope, intercept = np.polyfit(df['NDVI'].values,robust_y_convert,1)
                # renew the slope and intercept if the correlation coefficient is above 0.7
                if coef >= 0.7:
                    [renew_slope,renew_intercept] = [slope, intercept]
                else: pass
            else: 
                coef = -9999
                pass

            # gain index for soil and canopy pixel for each local domain
            index_veg = np.where(df['NDVI'] >= Veg_threshold)
            index_soil = np.where(df['NDVI'] <= Soil_threshold)
            # when the domain contains both vegetation and soil
            if len(index_veg[0]) > 0 and len(index_soil[0]) > 0:
                t_canopy[irow,icol] = df['Tr'][index_veg[0]].mean()
                t_soil[irow,icol] = df['Tr'][index_soil[0]].mean()
            # when the domain contains vegetation but no soil: estimate the soil temperature
            elif len(index_veg[0]) > 0 and len(index_soil[0]) == 0:
                t_canopy[irow,icol] = df['Tr'][index_veg[0]].mean()
                t_soil[irow,icol] = renew_slope*Soil_threshold + renew_intercept
            # when the domain contains soil but no vegetation: vegetation temperature is "NAN"
            elif len(index_veg[0]) == 0 and len(index_soil[0]) > 0:
                t_canopy[irow,icol] = np.nan
                t_soil[irow,icol] = df['Tr'][index_soil[0]].mean()
            # when the domain contains either pure soil or vegetation
            # estimate the soil and vegetation temperature
            elif len(index_veg[0]) == 0 and len(index_soil[0]) == 0:
                t_canopy[irow,icol] = renew_slope*Veg_threshold + renew_intercept
                t_soil[irow,icol] = renew_slope*Soil_threshold + renew_intercept
                if t_canopy[irow,icol] < -5:
                    t_canopy[irow,icol] = np.nan
                else: pass
                if t_soil[irow,icol] < -5:
                    t_soil[irow,icol] = np.nan
                else: pass
            t_coeff[irow,icol] = -1*coef

    tt_canopy = np.sqrt(np.sqrt(t_canopy.copy())) + 273.15
    tt_soil = np.sqrt(np.sqrt(t_soil.copy())) + 273.15
    tt_single_layer = (tt_canopy + tt_soil)/2
    # print("Canopy temperature is:",round(tt_canopy,2))
    # print("Soil temperature is:",round(tt_soil,2))

    # Write the separated temperature file
    driver = gdal.GetDriverByName('GTiff')
    ds = driver.Create(dir_output+"\\"+output_name_multiple, dims_LAI[1], dims_LAI[0], 3, gdal.GDT_Float32)
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
        os.remove(dir_output+"\\"+name_clip_opt)
        os.remove(dir_output+"\\"+name_clip_tr)
        os.remove(dir_output+"\\"+name_clip_dsm)
        os.remove(dir_output+"\\"+name_hillshade)
        os.remove(dir_output+"\\"+name_aggregate)
        print("Done!!! Middle products are deleted.")
    else:
        print("Done!!! Middle products are saved.")
        pass
    return