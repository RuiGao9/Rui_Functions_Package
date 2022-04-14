import arcpy
import arcview
import gdal
import numpy as np
import pandas as pd
import os

    
def TellExtent(InputFile):
    '''
    parameter InputFile: input could be Shpfile or image file
    return: the coordinates of the extent, 4 values
    '''
    
    Grid_Describe = arcpy.Describe(InputFile)
    Grid_Extent = Grid_Describe.extent
    ExtentShpfile = "{} {} {} {}".format(Grid_Extent.XMin, Grid_Extent.YMin, Grid_Extent.XMax, Grid_Extent.YMax)
    print("The extent is: " + ExtentShpfile)
    return(ExtentShpfile)
    
def TellResolution(Raster_File):
    '''
    To know the resolution of the raster file.
    :param Raster_File: put the raster file inside
    :return: the resolution (x,y) will be shown
    '''
    
    rastersize_x = arcpy.GetRasterProperties_management(Raster_File, "CELLSIZEX")
    rastersize_y = arcpy.GetRasterProperties_management(Raster_File, "CELLSIZEY")
    rastersize_x = round(float(rastersize_x.getOutput(0)), 3)
    rastersize_y = round(float(rastersize_y.getOutput(0)), 3)
    print("The horizontal resolution is " + str(rastersize_x) + "m and the vertical resolution is " + str(
        rastersize_y) + "m.")
    return (rastersize_x, rastersize_y)
    
def TellTheGeoInfo(InputRaster):
    '''
    parameter InputRaster: input image
    return: 1) dimension of the image, 2) geographic information, 3) projection information.
    '''
    
    fid=gdal.Open(InputRaster)
    input_img=fid.GetRasterBand(1).ReadAsArray()
    dims=input_img.shape
    print("Dimension of the data is:",dims[0],dims[1])
    # Read the GDAL GeoTransform to get the pixel size
    img_geo=fid.GetGeoTransform()
    img_prj=fid.GetProjection()
    fid=None
    return(dims,img_geo,img_prj)
    
def WriteTiffData(floder_out, name_out_file, ysize, xsize, Array_Content, geotransform, projection):
    '''
    Writing the array with geographic information.
    :param floder_out: directory that you wish to save this output
    :param name_out_file: a name for the output, like "Result"
    :param ysize: dimension on y
    :param xsize: dimension on x
    :param Array_Content: the array that you wish to save in the Tiff file
    :param geotransform: geotransform information
    :param projection: projection information
    :return: no return on screen, but an output in the folder
    '''
    import gdal
    
    driver = gdal.GetDriverByName('GTiff')
    print(floder_out + "\\" + name_out_file)
    new_tiff = driver.Create(floder_out + "\\" + name_out_file, xsize, ysize, 1, gdal.GDT_Float32)
    new_tiff.SetGeoTransform(geotransform)
    new_tiff.SetProjection(projection)
    new_tiff.GetRasterBand(1).WriteArray(Array_Content)
    new_tiff.FlushCache()  # Saves to disk
    new_tiff = None  # closes the file
    print("Done!!! Tiff data has been written.")
    return ()
    
def WriteTiffData_SingleOutput(floder_and_filename, ysize, xsize, Array_Content, geotransform, projection):
    '''
    Writing the array with geographic information.
    :param floder_out: directory that you wish to save this output
    :param name_out_file: a name for the output, like "Result"
    :param ysize: dimension on y
    :param xsize: dimension on x
    :param Array_Content: the array that you wish to save in the Tiff file
    :param geotransform: geotransform information
    :param projection: projection information
    :return: no return on screen, but an output in the folder
    '''
    
    driver = gdal.GetDriverByName('GTiff')
    print(floder_and_filename)
    new_tiff = driver.Create(floder_and_filename, xsize, ysize, 1, gdal.GDT_Float32)
    new_tiff.SetGeoTransform(geotransform)
    new_tiff.SetProjection(projection)
    new_tiff.GetRasterBand(1).WriteArray(Array_Content)
    new_tiff.FlushCache()  # Saves to disk
    new_tiff = None  # closes the file
    print("Done!!! Tiff data has been written.")
    return ()
    
def LST_CWSI(input_temp, input_boundary, output_dir, 
             output_temp_name, output_cwsi_name, NoDataValue):
    '''
    param input_temp: 2-layer temperature image: 1st layer is canopy temperature and 2nd layer is soil temperature.
    param input_boundary: research boundary. 
    param output_dir: a folder where you want to save the results.
    param output_temp_name: a name (with ".tif") for the temperature image within the research boundary.
    param output_cwsi_name: a name (with ".tif") for the CWSI image.
    return: CWSI map will be returned.
    '''
    
    # Extract raster only focus on the interesting area
    arcpy.gp.ExtractByMask_sa(input_temp, input_boundary, output_dir+"\\"+output_temp_name)
    # read the information from the image
    [res_tr_x,res_tr_y] = TellResolution(output_dir+"\\"+output_temp_name)
    extent = TellExtent(output_dir+"\\"+output_temp_name)
    [dims,img_geo,img_prj] = TellTheGeoInfo(output_dir+"\\"+output_temp_name)
    print("Image resolution is:",round(res_tr_x,2),round(res_tr_y,2))
    # print("Image geographic info:",img_geo)
    # print("Image projection info:",img_prj)
    array_temp = arcpy.RasterToNumPyArray(output_dir+"\\"+output_temp_name, nodata_to_value = NoDataValue)
    # calculate the CWSI
    # remember, the 1st layer is the canopy temperature
    [temp_max, temp_min] = [np.nanmax(array_temp[0,:,:]),np.nanmin(array_temp[0,:,:])]
    print("The maximum temperature within the image is",round(temp_max,2),"and the minimum is",round(temp_min,2))
    img_cwsi = (array_temp[0,:,:]-temp_min)/(temp_max-temp_min)
    # output the CWSI image
    WriteTiffData(output_dir, output_cwsi_name, dims[0], dims[1], img_cwsi, img_geo, img_prj)
    
    return print("Canopy temperature based CWSI image is generated!!!")
    
def Footprint_Digital_Results(footprint, tseb_r_1, tseb_r_2, temp_image, dir_out, 
                              lai_image, fc_image,
                              n_rn, n_h, n_le, n_g, n_t_et, pixel_size,
                              upper_boundary, lower_boundary, 
                              delete_tmp_files="Yes", single_layer_temp="Yes"):
    '''
    parameters:
    footprint: directory of the footprint image.
    tseb_r_1: directory of the TSEB result: multiple-layer image.
    tseb_r_2: directory of the TSEB ancillary result: multiple-layer image.
    temp_image: directory of the temperature image. it could be single-layer or multiple-layer image, but the "single_layer_temp" need to be set    correspondingly.
    dir_out: directory of the outputs from this scripts. they are transform results and they can be deleted.
    lai_image: directory of the LAI image.
    fc_image: directory of the fractional cover image.
    n_rn, n_h, n_le, n_g: the layer number of the net radiation, sensible heat flux, latent heat flux, and soil surface heat flux.
    n_t_et: the layer number of the ratio between canopy latent heat flux and total latent heat flux.
    pixel_size: the pixel size (e.g., 3.6 meter by 3.6 meter).
    upper_boundary: the upper threshold for all fluxes at one pixel which does not make sense, e.g., 10,000 W/m2 for LE at one pixel.
    lower_boundary: the lower threshold for all fluxes at one pixel which does not make sense, e.g., -1,500 W/m2 for LE at one pixel.
    delete_tmp_files: Default is "Yes", and this means the temporary (middle products) files will be deleted at the end. Any other input (string) results in saving the temporary files.
    single_layer_temp: Default is "Yes", and this means a single layer temperature image. Other input (string) must be multiple layers temperature image, 
                        and the 1st and 2nd layer must be the canopy and soil temperature, respectively.

    return:
    Net radiation, sensible heat flux, latent heat flux, soil surface heat flux, canopy latent heat flux, 
            LAI, single-layer mean temperature, canopy temperature, and soil temperature within the footprint area.
    '''
    
    cellsize = str(pixel_size)+" "+str(pixel_size)
    arcpy.Resample_management(in_raster=footprint, 
                              out_raster=dir_out+"\\footprint_resample.tif", 
                              cell_size=cellsize, 
                              resampling_type="CUBIC")

    Grid_Describe = arcpy.Describe(tseb_r_1)
    Grid_Extent = Grid_Describe.extent
    extent = "{} {} {} {}".format(Grid_Extent.XMin, Grid_Extent.YMin, Grid_Extent.XMax, Grid_Extent.YMax)

    arcpy.Clip_management(in_raster=dir_out+"\\footprint_resample.tif", 
                          rectangle=extent, 
                          out_raster=dir_out+"\\footprint_clip.tif", 
                          in_template_dataset=tseb_r_1, 
                          nodata_value="0.000000e+00", 
                          clipping_geometry="NONE", maintain_clipping_extent="MAINTAIN_EXTENT")

    raster_footprint = arcpy.RasterToNumPyArray(dir_out+"\\footprint_clip.tif", nodata_to_value=-9999)
    raster_footprint[raster_footprint>1] = np.nan
    raster_footprint[raster_footprint<0] = np.nan
    raster_tseb = arcpy.RasterToNumPyArray(tseb_r_1, nodata_to_value=-9999)
    raster_tseb_ancillary = arcpy.RasterToNumPyArray(tseb_r_2, nodata_to_value=np.nan)
    raster_lai = arcpy.RasterToNumPyArray(lai_image, nodata_to_value=np.nan)
    raster_fc = arcpy.RasterToNumPyArray(fc_image, nodata_to_value=np.nan)

    raster_rn = raster_tseb[n_rn,:,:]
    raster_rn[raster_rn>upper_boundary] = np.nan
    raster_rn[raster_rn<lower_boundary] = np.nan
    out_rn = raster_rn*raster_footprint

    raster_h = raster_tseb[n_h,:,:]
    raster_h[raster_h>upper_boundary] = np.nan
    raster_h[raster_h<lower_boundary] = np.nan
    out_h = raster_h*raster_footprint

    raster_le = raster_tseb[n_le,:,:]
    raster_le[raster_le>upper_boundary] = np.nan
    raster_le[raster_le<lower_boundary] = np.nan
    out_le = raster_le*raster_footprint

    raster_g = raster_tseb[n_g,:,:]
    raster_g[raster_g>upper_boundary] = np.nan
    raster_g[raster_g<lower_boundary] = np.nan
    out_g = raster_g*raster_footprint

    raster_t = raster_tseb_ancillary[n_t_et,:,:]
    raster_t[raster_t>1] = np.nan
    raster_t[raster_t<0] = np.nan
    out_t = out_le*raster_t

    raster_tet = raster_footprint*0+1
    out_tet = raster_tet*raster_tseb_ancillary[n_t_et,:,:]
    out_tet = np.nanmean(out_tet)    

    raster_lai[raster_lai>5] = np.nan
    raster_lai[raster_lai<0] = np.nan
    out_lai = raster_lai*(raster_footprint*0+1)

    raster_fc[raster_fc>1] = np.nan
    raster_fc[raster_fc<0] = np.nan
    out_fc = raster_fc*(raster_footprint*0+1)

    if single_layer_temp == "Yes":
        raster_temp = arcpy.RasterToNumPyArray(temp_image, nodata_to_value=np.nan)
        raster_temp[raster_temp>350] = np.nan # 76.85 C
        raster_temp[raster_temp<270] = np.nan # -3.15 C
        out_temp = raster_temp*(raster_footprint*0+1)

        out_rn = np.nansum(out_rn)
        out_h = np.nansum(out_h)
        out_le = np.nansum(out_le)
        out_g = np.nansum(out_g)
        out_t = np.nansum(out_t)
        out_lai = np.nanmean(out_lai)
        out_fc = np.nanmean(out_fc)
        out_temp = np.nanmean(out_temp)
        out_temp_canopy = np.nan
        out_temp_soil = np.nan
        print("Rn - Net radiation:",round(out_rn,3))
        print("H - Sensible heat flux:",round(out_h,3))
        print("LE - Latent heat flux:",round(out_le,3))
        print("G - Soil surface heat flux:",round(out_g,3))
        print("T - Canopy latent heat flux:",round(out_t,3))
        print("ET partitioning:",round(out_tet,3))
        print("\nLAI:",round(out_lai,3))
        print("Fractional cover:",round(out_fc,3))
        print("Temperautre:",round(out_temp,3),"K")
    else: 
        # canopy, soil temperature are required at the 1st and 2nd layer of the image
        raster_temp = arcpy.RasterToNumPyArray(temp_image, nodata_to_value=np.nan)
        raster_temp_canopy = raster_temp[0,:,:]
        raster_temp_soil = raster_temp[1,:,:]

        raster_temp_canopy[raster_temp_canopy>350] = np.nan # 76.85 C
        raster_temp_canopy[raster_temp_canopy<270] = np.nan # -3.15 C
        out_temp_canopy = raster_temp_canopy*(raster_footprint*0+1)

        raster_temp_soil[raster_temp_soil>350] = np.nan # 76.85 C
        raster_temp_soil[raster_temp_soil<270] = np.nan # -3.15 C
        out_temp_soil = raster_temp_soil*(raster_footprint*0+1)

        out_rn = np.nansum(out_rn)
        out_h = np.nansum(out_h)
        out_le = np.nansum(out_le)
        out_g = np.nansum(out_g)
        out_t = np.nansum(out_t)
        out_lai = np.nanmean(out_lai)
        out_fc = np.nanmean(out_fc)
        out_temp_canopy = np.nanmean(out_temp_canopy)
        out_temp_soil = np.nanmean(out_temp_soil)
        out_temp = np.nan
        print("Rn - Net radiation:",round(out_rn,3))
        print("H - Sensible heat flux:",round(out_h,3))
        print("LE - Latent heat flux:",round(out_le,3))
        print("G - Soil surface heat flux:",round(out_g,3))
        print("T - Canopy latent heat flux:",round(out_t,3))
        print("ET partitioning:",round(out_tet,3))
        print("\nLAI:",round(out_lai,3))
        print("Fractional cover:",round(out_fc,3))
        print("Canopy temperautre:",round(out_temp_canopy,3),"K")
        print("Soil temperautre:",round(out_temp_soil,3),"K")

    if delete_tmp_files == "Yes":
        os.remove(dir_out+"\\footprint_resample.tif")
        os.remove(dir_out+"\\footprint_clip.tif")
    else:
        print("Temporary files are saved in the output folder.")

    return(out_rn, out_h, out_le, out_g, out_t, out_tet, out_lai, out_fc, out_temp, out_temp_canopy, out_temp_soil)
    
