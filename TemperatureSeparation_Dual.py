def TemperatureSeparation_Dual(dir_LAI, dir_RGBNIR, dir_Tr, dir_DSM, 
                               NoDataValue, Veg_threshold, Soil_threshold,
                               band_R, band_NIR, 
                               cellsize_resample, cellsize_grid,
                               dir_output, output_name_multiple, output_name_single,
                               Azimuth, Altitude):
    # --- --- --- ---
    # Single layer temperature processing
    # 1. Upscale the temperature image
    name_cal_tr = "tr_calculator.tif"
    Array_Tr_upscale = arcpy.RasterToNumPyArray(dir_Tr, nodata_to_value=NoDataValue)
    [res_tr_x,res_tr_y] = TellResolution(dir_Tr)
    ## Get the information from temperature map for data output
    fid_tr_upscale=gdal.Open(dir_Tr)
    input_tr_upscale=fid_tr_upscale.GetRasterBand(1).ReadAsArray()
    dims_tr_upscale=input_tr_upscale.shape
    tr_geo_upscale=fid_tr_upscale.GetGeoTransform()
    tr_prj_upscale=fid_tr_upscale.GetProjection()
    fid_tr_upscale=None
    ## Stefan-Boltzmann Law
    Array_Tr_upscale = (Array_Tr_upscale+273.15) ** 4
    dims_Tr = Array_Tr_upscale.shape
    print("Column of the temperature image:",dims_tr_upscale[0],"Row of the temperature image:",dims_tr_upscale[1])
    print("Resolution of the pixel is",round(res_tr_x,2),"and",round(res_tr_y,2))
    WriteTiffData(dir_output, name_cal_tr, dims_tr_upscale[0], dims_tr_upscale[1], Array_Tr_upscale, tr_geo_upscale, tr_prj_upscale)

    # 2. Resample the temperature image from 0.6 m to 3.6 m
    name_resample_tr = "tr_resample.tif"
    tmp_tr = arcpy.ia.Resample(dir_output + "\\" + name_cal_tr, "Average", res_tr_x, cellsize_grid)
    tmp_tr.save(dir_output+"\\"+name_resample_tr)

    # 3. Downscale the temperature image
    name_cal_downscale_tr = "tr_calculator_downscale.tif"
    Array_Tr_downscale = arcpy.RasterToNumPyArray(dir_output+"\\"+name_resample_tr, nodata_to_value=NoDataValue)
    Array_Tr_downscale = np.sqrt(np.sqrt(Array_Tr_downscale))
    [res_tr_downscale_x,res_tr_downscale_y] = TellResolution(dir_output+"\\"+name_resample_tr)
    ## Get the information from temperature map for data output
    fid_tr_downscale=gdal.Open(dir_output+"\\"+name_resample_tr)
    input_tr_downscale=fid_tr_downscale.GetRasterBand(1).ReadAsArray()
    dims_tr_downscale=input_tr_downscale.shape
    tr_downscale_geo=fid_tr_downscale.GetGeoTransform()
    tr_downscale_prj=fid_tr_downscale.GetProjection()
    fid_tr_downscale=None
    print("Column of the temperature image:",dims_tr_downscale[0],"Row of the temperature image:",dims_tr_downscale[1])
    print("Resolution of the pixel is",round(res_tr_downscale_x,2),"and",round(res_tr_downscale_y,2))
    WriteTiffData(dir_output, name_cal_downscale_tr, dims_tr_downscale[0], dims_tr_downscale[1], Array_Tr_downscale, tr_downscale_geo, tr_downscale_prj)

    # 4. Clip the raster as LAI
    # Clip the temperature image and the clipped one is the final result
    extent = TellExtent(dir_LAI)
    arcpy.Clip_management(in_raster=dir_output+"\\"+name_cal_downscale_tr, 
                          rectangle=extent, 
                          out_raster=dir_output+"\\"+output_name_single, 
                          in_template_dataset=dir_LAI, 
                          nodata_value=NoDataValue, 
                          clipping_geometry="NONE", 
                          maintain_clipping_extent="MAINTAIN_EXTENT")  

    # --- --- --- ---
    # Muitlple layer temperature processing
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
    # print(Array_LAI.shape, Array_RGBNIR.shape, Array_R.shape, Array_NIR.shape, Array_Tr.shape)
    # print(round(np.nanmax(Array_NDVI),1),round(np.nanmin(Array_NDVI),1))

    # Stefan-Boltzmann Law
    Array_Tr = (Array_Tr+273.15) ** 4

    dims_LAI = Array_LAI.shape
    print("Column of the LAI:",dims_LAI[0],"\nRow of the LAI:",dims_LAI[1])
    dims_NDVI = Array_NDVI.shape
    print("Column of the spectral data:",dims_NDVI[0],"\nRow of the spectral data:",dims_NDVI[1])
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
    # Initial arrays
    t_canopy = np.empty((dims_LAI[0],dims_LAI[1]))
    t_canopy[:] = np.nan
    t_soil = np.empty((dims_LAI[0],dims_LAI[1]))
    t_soil[:] = np.nan
    t_coeff = np.empty((dims_LAI[0],dims_LAI[1]))
    t_coeff[:] = np.nan
    print("Dimension of the canopy temperature is:",t_canopy.shape[0],t_canopy.shape[1])
    print("Dimension of the soil temperature is:",t_soil.shape[0],t_soil.shape[1])

    # Initial parameters which is convenient for model checking
    [renew_slope,renew_intercept,
     model_score_1,model_score_2] = [NoDataValue,NoDataValue,[],[]]
    # create an initial robust model
    data = [[0.1,0],[0.5,0],[1,0],[10,0],[100,0]]
    df_temp = pd.DataFrame(data, columns = ['X','Y'])
    model_final = linear_model.RANSACRegressor()
    model_final.fit(df_temp[['X']],df_temp[['Y']])
    reg_final = model_final.fit(df_temp[['X']],df_temp[['Y']])
    model_temp_score = round(reg_final.score(df_temp[['X']],df_temp[['Y']]),2)
    print("The score of the initial model is",model_temp_score)

    for irow in range(dims_LAI[0]):
        start_row = irow * hor_pixel
        end_row = start_row + (hor_pixel)
        for icol in range(dims_LAI[1]):
            start_col = icol * ver_pixel
            end_col = start_col + (ver_pixel)
#             print("\n\nRow number is:",irow,"Column number is:",icol)

            # Step 1. Using the hillshade to eliminate the shadow pixel
            local_NDVI = Array_NDVI[start_row:end_row,start_col:end_col]
            local_NDVI[local_NDVI < 0] = NoDataValue
            local_Tr = Array_Tr[start_row:end_row,start_col:end_col] # temperature power 4
            local_Shadow = Array_Shadow[start_row:end_row,start_col:end_col]
            num_zero = np.count_nonzero(local_Shadow==0)
            local_NDVI = local_NDVI*local_Shadow
            tmp_NDVI = local_NDVI.reshape(-1)
            tmp_NDVI[tmp_NDVI<=0] = np.nan
            tmp_Tr = local_Tr.reshape(-1)
            tmp_Tr = np.sqrt(np.sqrt(tmp_Tr)) # unit in degree K
            # creat a data frame without the impact of shadow pixels
            df = pd.DataFrame()
            df = pd.DataFrame({'NDVI': tmp_NDVI,'Tr': tmp_Tr})
            # delete infinite and nan values
#             df_soil = df_copy[(df_copy["NDVI"] <= Soil_threshold)]
            df[df['NDVI']<0] = np.nan
            df[df['NDVI']>1] = np.nan
            # The temperature within [273.15,350] is recognized as valid values
            df[df['Tr']<273.15] = np.nan
            df[df['Tr']>350] = np.nan
            df = df.dropna()
            df = df.drop_duplicates()
            df = df.apply(pd.to_numeric, errors='coerce')

            # Step 2. 1st round seeing the relationship
            # 25% of the pixels in the domain is set as a threshold, and this is based on personal empirical experience
            if len(df) > 9:
                # Robust linear relationship
                robust_model_1 = linear_model.RANSACRegressor()
                robust_model_1.fit(df[['NDVI']],df[['Tr']])
                reg = robust_model_1.fit(df[['NDVI']],df[['Tr']])
                score_1 = round(reg.score(df[['NDVI']],df[['Tr']]),2)
#                 print("Score of the 1st model is:",score_1)
#                 print("Estimated soil temperature is:",robust_model_1.predict([[Soil_threshold]]),
#                       "\nEstimated canopy temperature is:",robust_model_1.predict([[Veg_threshold]]))
                model_score_1.append(score_1)

                # To get the coefficient for the 3rd layer output information
                robust_y = robust_model_1.predict(df[['NDVI']])                        
                # processing for robust_y and using a general algorithm to describe the robust linear relationship
                [num, robust_y_convert] = [len(robust_y),[]]
                for irou in range(0,num):
                    robust_y_convert.append(robust_y[irou][0])
                # correlation coefficient between the temperature predictions and observations
                coef_1 = gfit(df['Tr'].values,np.asarray(robust_y_convert),type_statistic='4',residual='No')
#                 print("Valid pixel number in the 1st round is:",len(df['NDVI']))
#                 print("Correlation coefficient between prediction and observation is:",round(coef_1,2))
                # getting the slope and intercept from the robust linear regression
                slope, intercept = np.polyfit(df['NDVI'].values,robust_y_convert,1)
#                 print("Slope and intercept calculated based on the robust linear regression model are:\n",round(slope,2),"and",round(intercept,2))
            else: 
                coef = 99999
                pass

            # Step 3: 2nd round seeing the relationship
            # 25% of the pixels within the domain is set as a threshold, and this is based on personal empirical experience
            if len(df) > 9:
#                 # Plotting a scatter plot to describe the 1st robust linear relationship
#                 plt.plot(df['NDVI'],df['Tr'],'ro')
#                 plt.plot(df['NDVI'],robust_y)
#                 plt.xlabel('NDVI')
#                 plt.ylabel('Temperature K')
#                 plt.show()
                # Check the temperature distribution for veg and soil
                df_soil = df[(df["NDVI"] <= Soil_threshold)]
                df_veg = df[(df["NDVI"] >= Veg_threshold)]
                df_mid = df[(df["NDVI"] > Soil_threshold) &
                            (df["NDVI"] < Veg_threshold)]
                # Experiment which is not a part of this program - Sorting NDVI and Tr (NDVI is increasing and Tr is decreasing)
                df_mid_sorted = pd.DataFrame({'NDVI': df_mid['NDVI'].sort_values(ascending=True).values,'Tr': df_mid['Tr'].sort_values(ascending=False).values})
#                 # Box plot checking those three categories
#                 fig, ([ax1, ax2, ax3]) = plt.subplots(1,3)
#                 ax1.boxplot(df_soil["Tr"])
#                 ax2.boxplot(df_mid["Tr"])
#                 ax3.boxplot(df_veg["Tr"])
#                 plt.show()
#                 print("The number of soil pixel (<=0.40):",len(df_soil['Tr']))
#                 print("The number of veg pixel  (>=0.70):",len(df_veg['Tr']))

                # Elimiting pixels to gain the NDVI-Tr relationship
                # Rules:
                # 1. Ignore the vegetation point whose temperature is above the 50 percentile of the original vegetation points
                # Note, initially, we decide to ignore the temperature which is above 75 percentile, 
                # but 50 percentile seems better in the later process.
                temp_veg_table = df_veg.describe()
#                 print("75 percentile veg temperature:",round(temp_veg_table.values[6,1],2))
#                 print("50 percentile veg temperature:",round(temp_veg_table.values[5,1],2))
#                 print("25 percentile veg temperature:",round(temp_veg_table.values[4,1],2))
                veg_75 = temp_veg_table.values[6,1]
                veg_50 = temp_veg_table.values[5,1]
                veg_25 = temp_veg_table.values[4,1]
                # drop the temperature which is above the 50 percentile in the vegetation area
                df_veg[df_veg["Tr"] >= veg_50] = np.nan
                df_veg = df_veg.dropna()
                df_veg = df_veg.apply(pd.to_numeric, errors='coerce')
                df_soil = df_soil.dropna()
                df_soil = df_soil.apply(pd.to_numeric, errors='coerce')
#                 print("The number of shrinked vegetation pixels is:",len(df_veg["Tr"]))
#                 print("The number of shrinked soil pixels is:",len(df_soil["Tr"]))
                # concat those three pieces and then gain the robust linear relationship again
                df = pd.concat([df_soil,df_mid,df_veg])

                # Put this condition [len(df)>9] here in case the number of valid records is not sufficient
                # 25% of the pixel is set as a threshold, and this is based on personal empirical experience
                if len(df) > 9:
                    robust_model_2 = linear_model.RANSACRegressor()
                    robust_model_2.fit(df[['NDVI']],df[['Tr']])
                    reg_2 = robust_model_2.fit(df[['NDVI']],df[['Tr']])
                    score_2 = round(reg_2.score(df[['NDVI']],df[['Tr']]),2)
#                     print("\nScore of the 2nd model is:",score_2)
#                     print("Estimated soil temperature is:",robust_model_2.predict([[Soil_threshold]]),
#                           "\nEstimated canopy temperature is:",robust_model_2.predict([[Veg_threshold]]))
                    model_score_2.append(score_2)
                    # To get the coefficient for the 3rd layer output information
                    robust_y_2 = robust_model_2.predict(df[['NDVI']]) 
                    # processing for robust_y_2 and using a general algorithm to describe the robust linear relationship
                    [num_2, robust_y_convert_2] = [len(robust_y_2),[]]
                    for irou_2 in range(0,num_2):
                        robust_y_convert_2.append(robust_y_2[irou_2][0])
                    # correlation coefficient between the temperature predictions and observations
                    coef_2 = gfit(df['Tr'].values,np.asarray(robust_y_convert_2),type_statistic='4',residual='No')
#                     print("Valid pixel number in the 2nd round is:",len(df['NDVI']))
#                     print("Correlation coefficient between prediction and observation is:",round(coef_2,2))
                    # getting the slope and intercept from the robust linear regression
                    slope_2, intercept_2 = np.polyfit(df['NDVI'].values,robust_y_convert_2,1)
#                     print("Slope and intercept calculated based on the robust linear regression model are:\n",round(slope_2,2),"and",round(intercept_2,2))
                    # remember the parameters based on the better robust relationship
                    # 0.5 is an empirical parameter based on the interpretation of the data
                    if score_2 > score_1 and score_2 >= 0.5 and len(df_veg) >= 1 and len(df_mid) >= 1:
                        # remember the 2nd robust model
                        model_final = robust_model_2
                        coef = coef_2
#                         print("The 2nd robust model is remembered instantaneously")
                    elif score_1 > score_2 and score_1 >= 0.5 and len(df_veg) >= 1 and len(df_mid) >= 1:
                        # remember the 1st robust model
                        model_final = robust_model_1
                        coef = coef_1
#                         print("The 1st robust model is remembered instantaneously")
                    else: 
                        coef = 99999
                        pass
#                     # plot scatter plot to describe the 2nd robust linear relationship
#                     plt.plot(df['NDVI'],df['Tr'],'ro')
#                     plt.plot(df['NDVI'],robust_y_2)
#                     plt.xlabel('NDVI')
#                     plt.ylabel('Temperature K')
#                     plt.show()
                else:
                    t_canopy[irow,icol] = model_final.predict([[Veg_threshold]])
                    t_soil[irow,icol] = model_final.predict([[Soil_threshold]])
#                     print(type(t_canopy[irow,icol]))
        #             print("Valid pixels are not sufficient for the algorithm, and estimations are made.")
        #             print("Canopy temperature is estimated:",round(t_canopy[irow,icol],2))
        #             print("Soil temperature is estimated:",round(t_soil[irow,icol],2))
        #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2)):
            else:
                t_canopy[irow,icol] = model_final.predict([[Veg_threshold]])
                t_soil[irow,icol] = model_final.predict([[Soil_threshold]])
    #             print("Valid pixels are not sufficient for the algorithm, and estimations are made.")
    #             print("Canopy temperature is estimated:",round(t_canopy[irow,icol],2))
    #             print("Soil temperature is estimated:",round(t_soil[irow,icol],2))
    #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2)):

            # canopy and soil temperature arrangement based on the above result
            df_copy = df.copy(deep=True)
            df_soil = df_copy[(df_copy["NDVI"] <= Soil_threshold)]
            df_veg = df_copy[(df_copy["NDVI"] >= Veg_threshold)]
            df_mid = df_copy[(df_copy["NDVI"] > Soil_threshold) &
                             (df_copy["NDVI"] < Veg_threshold)]
            [pixel_canopy,pixel_soil] = [len(df_veg),len(df_soil)]
            # when the domain contains both vegetation and soil
            if pixel_canopy > 0 and pixel_soil > 0:
                t_canopy[irow,icol] = df_veg['Tr'].mean()
                t_soil[irow,icol] = df_soil['Tr'].mean()
#                 print("Canopy temperature is gained directly:",round(t_canopy[irow,icol],2))
#                 print("Soil temperature is gained directly:",round(t_soil[irow,icol],2))
    #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2))
            # when the domain contains vegetation but no soil: estimate the soil temperature
            elif pixel_canopy > 0 and pixel_soil == 0:
                t_canopy[irow,icol] = df_veg['Tr'].mean()
                t_soil[irow,icol] = model_final.predict([[Soil_threshold]])
#                 print("Canopy temperature is gained directly:",round(t_canopy[irow,icol],2))
#                 print("Soil temperature is estimated:",round(t_soil[irow,icol],2))
    #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2))
            # when the domain contains soil but no vegetation: vegetation temperature is "NAN"
            elif pixel_canopy == 0 and pixel_soil > 0:
                t_canopy[irow,icol] = np.nan
                t_soil[irow,icol] = df_soil['Tr'].mean()
#                 print("Canopy temperature is estimated:",round(t_canopy[irow,icol],2))
#                 print("Soil temperature is gained directly:",round(t_soil[irow,icol],2))
    #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2))
            # when the domain contains either pure soil or vegetation
            # estimate the soil and vegetation temperature
            elif pixel_canopy == 0 and pixel_soil == 0:
                t_canopy[irow,icol] = model_final.predict([[Veg_threshold]])
                t_soil[irow,icol] = model_final.predict([[Soil_threshold]])
#                 print("Canopy temperature is estimated:",round(t_canopy[irow,icol],2))
#                 print("Soil temperature is estimated:",round(t_soil[irow,icol],2))
    #             print("The used slope and intercept for this round is:",round(renew_slope,2),round(renew_intercept,2))
            t_coeff[irow,icol] = -1*coef
            # debugging
#             diff_tmp = t_canopy[irow,icol] - t_soil[irow,icol]
#             if diff_tmp > 0:
#                 print("The canopy temperature is above the soil temperature!!!")
#             else: pass
    tt_canopy = t_canopy.copy()
    tt_soil = t_soil.copy()

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
    
    return()