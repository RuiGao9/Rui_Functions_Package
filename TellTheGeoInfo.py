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