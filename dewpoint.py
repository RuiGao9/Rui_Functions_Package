def dewpoint(rh,temp):
    '''
    This function is built for dew point calculation.
    
    Parameters:
    rh: relative humidity, from sensor (%)
    temp: corresponding temperature

    Return:
    Dew point in degree C will be returned.
    '''
	
    import math
    a = 17.27
    b = 237.7
    alpha = a*temp/(b+temp)+math.log(rh/100)
    td = b*alpha/(a-alpha)
    
    return(td)