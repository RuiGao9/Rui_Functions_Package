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
    import arcpy

    driver = gdal.GetDriverByName('GTiff')
    print(floder_out + "\\" + name_out_file + ".tif")
    new_tiff = driver.Create(floder_out + "\\" + name_out_file + ".tif", xsize, ysize, 1, gdal.GDT_Float32)
    new_tiff.SetGeoTransform(geotransform)
    new_tiff.SetProjection(projection)
    new_tiff.GetRasterBand(1).WriteArray(Array_Content)
    new_tiff.FlushCache()  # Saves to disk
    new_tiff = None  # closes the file
    print("Done!!! Tiff data has been written.")
    return ()