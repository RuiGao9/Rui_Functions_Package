class ArcGIS_Rui:
    def __init__(self,Raster_File=""):
        self.Raster_File = Raster_File
        
    
    import arcpy
    import arcview
    import numpy as np
    
    def TellResolution(self,Raster_File=""):
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
    
    def TellExtent(self,Raster_File=""):
        '''
        parameter InputFile: input could be Shpfile or image file
        return: the coordinates of the extent, 4 values
        '''
        Grid_Describe = arcpy.Describe(Raster_File)
        Grid_Extent = Grid_Describe.extent
        ExtentShpfile = "{} {} {} {}".format(Grid_Extent.XMin, Grid_Extent.YMin, Grid_Extent.XMax, Grid_Extent.YMax)
        print("The extent is: " + ExtentShpfile)
        return(ExtentShpfile)
    
    def TellTheGeoInfo(self,Raster_File=""):
        '''
        parameter InputRaster: input image
        return: 1) dimension of the image, 2) geographic information, 3) projection information.
        '''
        fid=gdal.Open(Raster_File)
        input_img=fid.GetRasterBand(1).ReadAsArray()
        dims=input_img.shape
        print("Dimension of the data is:",dims[0],dims[1])
        # Read the GDAL GeoTransform to get the pixel size
        img_geo=fid.GetGeoTransform()
        img_prj=fid.GetProjection()
        fid=None
        return(dims,img_geo,img_prj)