import arcpy
import arcview
import gdal

    
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
             output_temp_name, output_cwsi_name):
    '''
    param input_temp: temperature image.
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
    array_temp = arcpy.RasterToNumPyArray(output_dir+"\\"+output_temp_name, nodata_to_value = NoDataValue)
    # calculate the CWSI
    [temp_max, temp_min] = [np.nanmax(array_temp),np.nanmin(array_temp)]
    print("The maximum temperature within the image is",round(temp_max,2),"and the minimum is",round(temp_min,2))
    img_cwsi = (array_temp-temp_min)/(temp_max-temp_min)
    # output the CWSI image
    WriteTiffData(output_dir, output_cwsi_name, dims[0], dims[1], img_cwsi, img_geo, img_prj)
    
    return print("CWSI is written!!!")