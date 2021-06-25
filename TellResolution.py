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