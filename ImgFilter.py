def ImgFilter(dir_img,
              name_out_file,dir_output,
              value_min,value_max,value_nan):
    '''
    Eliminate the "invalid" pixel value.
    :param dir_img: the image you want to process.
    :param name_out_file: a name with ".tif" (the extension) for the output.
    :param dir_output: a folder where you want to save the output.
    :param value_min: the lowest boundary you defined for the image - any value below this boundary will be recognized as "value_nan".
    :param value_max: the highest boundary you defined for the image.
    :param value_nan: a value to define those "invalid" values.
    :return: a new image without the pixel value beyond the boundary.
    '''
    import arcpy
    import numpy as np
    import sys
    
    def TellExtent(InputFile):
        '''
        parameter InputFile: input could be Shpfile or image file
        return: the coordinates of the extent, 4 values
        '''
        import arcpy
        import arcview
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
        import arcpy
        import arcview

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
        import arcpy
        import gdal
        
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
        import arcpy
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
    
    img = arcpy.RasterToNumPyArray(dir_img, nodata_to_value=value_nan)
    extent = TellExtent(dir_img)
    [res_x,res_y] = TellResolution(dir_img)
    [dims,img_geo,img_prj] = TellTheGeoInfo(dir_img)
    print("The minimum value is",round(img.min(),2),"and maximum value is",round(img.max(),2))
    img[img<value_min] = value_nan
    img[img>value_max] = value_nan
    print("After processing, \nthe minimum value is", round(img.min(),2), "and the maximum value is", round(img.max(),2))
    WriteTiffData(dir_output,name_out_file,dims[0],dims[1],img,img_geo,img_prj)
    
    return()