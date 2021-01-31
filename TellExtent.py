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