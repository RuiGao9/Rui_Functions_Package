def VegetationIndex(OpticalImage,
                    band_R,band_G,band_B,band_NIR,
                    novalue,VegetationIndex="NDVI"):
    
    import arcpy
    import numpy as np
    
    raster_rgbnir = arcpy.RasterToNumPyArray(OpticalImage, nodata_to_value=novalue)
    raster_r = raster_rgbnir[band_R,:,:]
    raster_g = raster_rgbnir[band_G,:,:]
    raster_b = raster_rgbnir[band_B,:,:]
    raster_nir = raster_rgbnir[band_NIR,:,:]
    
    if VegetationIndex == "NDVI":
        VI = (raster_nir-raster_r)/(raster_nir+raster_r)
    elif VegetationIndex == "SR":
        VI = raster_nir/raster_r
    elif VegetationIndex == "NDWI":
        VI = (raster_g-raster_nir)/(raster_g+raster_nir)
    elif VegetationIndex == "GNDVI":
        VI = (raster_nir-raster_g)/(raster_g+raster_nir)
    elif VegetationIndex == "SAVI":
        VI = (raster_nir-raster_r)*1.5/(raster_r+raster_nir+0.5)
    elif VegetationIndex == "PVI":
        VI = (raster_nir-0.3*raster_r-0.5)/(np.sqrt(1+0.009))
    elif VegetationIndex == "EVI":
        VI = 2.5*(raster_nir-raster_r)/(raster_nir+6*raster_r-7.5*raster_b+1)
    elif VegetationIndex == "CIg":
        VI = raster_nir/raster_g-1
    elif VegetationIndex == "MTVI2":
        VI = 1.5*(1.2*(raster_nir-raster_g)-2.5*(raster_r-raster_g))/np.sqrt((2*raster_nir+1)**2-(-5*np.sqrt(raster_r))-0.5)
    elif VegetationIndex == "MSAVI":
        VI = 0.5*(2*(raster_nir+1)-np.sqrt((2*raster_nir+1)**2-8*(raster_nir-raster_r)))
    elif VegetationIndex == "TSAVI":
        VI = (0.33*(raster_nir-0.33*raster_r-0.5))/(0.5*raster_nir+raster_r-0.5*0.33+1.5*(1+0.33**2))
    elif VegetationIndex == "VARI":
        VI = (raster_g-raster_r)/(raster_g+raster_r-raster_b)
    elif VegetationIndex == "IronOxide":
        VI = raster_r/raster_b
    elif VegetationIndex == "BAI":
        VI = 1/((0.1-raster_r)**2+(0.06-raster_nir)**2)
    elif VegetationIndex == "GEMI":
        VI = (2*(raster_nir**2-raster_r**2)+1.5*raster_nir+0.5*raster_r)/(raster_nir+raster_r+0.5)
        VI = VI*(1-0.25*VI)-((raster_r-0.125)/(1-raster_r))
    
    return(VI)
