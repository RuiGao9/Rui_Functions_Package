def calculate_doy(year, month, day):
    '''
    Calculate DOY based on year, month, and day.
    
    '''
    from datetime import datetime
    
    date_object = datetime(year, month, day)
    doy = date_object.timetuple().tm_yday

    return doy