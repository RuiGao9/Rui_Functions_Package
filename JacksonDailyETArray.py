def JacksonDailyETArray(LE,Latitude,DOY,t,Ratio_W_per_SquareMeter_2_mm_per_day):
    '''
    Jackson's method can refer to https://www.sciencedirect.com/science/article/pii/0378377483900951
    Goal: estimate daily ET based on Jackson's approach via putting LE (an array with unit in W/m2). The output is the estimated an ET array with unit in mm/day.
    parameter LE: it is an array, W/m2
    parameter Latitude: it is decimal
    parameter DOY: the day of the year
    parameter t: the time period between the sunrise and current time
    parameter Ratio_W_per_SquareMeter_2_mm_per_day: the value changing W/m2 to mm/day
    return ET_Daily: return an array, and the unit is mm/day
    '''
    import math
    import numpy as np
    
    pi = math.pi
    a = 12 - 5.69**(-12)*Latitude - 2.02*10**(-4)*Latitude**2 + 8.25*10**(-6)*Latitude**3 - 3.15*10**(-7)*Latitude**4
    b = 0.123*Latitude - 3.10*10**(-4)*Latitude**2 + 8*10**(-7)*Latitude**3 - 4.99*10**(-7)*Latitude**4
    N = 0.945 * (a + b * math.sin(pi * (DOY+10)/365)**2);
    J = 2*N/(pi*math.sin(pi*t/N))
    LE_mm_day = LE*Ratio_W_per_SquareMeter_2_mm_per_day
    ET_Daily = J*LE_mm_day
    print('Letent Heat flux is measured',str(round(t,2)),
          'hr after the sunrise.\nThe time period between sunrise and sunset is:',str(round(N,2)),
          'hr.\nParameter a is:',str(round(a,2)),
          '\nParameter b is:',str(round(b,2)),
          '\nCoefficient J is:',str(round(J,2)))
    return ET_Daily