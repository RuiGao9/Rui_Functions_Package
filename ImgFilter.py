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