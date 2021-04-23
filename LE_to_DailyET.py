def LE_to_DailyET(LE,Approach='Rs',
                  Rsi=800,Rsd=8000,
                  Rni=500,Rnd=8000,
                  Gi=500,Gd=4000,
                  c=3600,rou=1000,lamb=2259.36):

    '''
    parameter LE: the instantent latent heat flux in W/m^2.
    parameter Approach: this one is string and it is restricted in 'Rs', 'EF', and 'Rn-Rs'.
    
    parameter Rsi: the instantent solar radiation in W/m^2.
    parameter Rsd: the summation of the instantent solar radiation in W/m^2.
    parameter Rni: the instantent net radiation in W/m^2.
    parameter Rnd: the summation of the instantent net radiation in W/m^2.
    parameter Gi: the instantent soil heat flux in W/m^2.
    parameter Gd: the summation of the instantent soil heat flux in W/m^2.
    
    parameter c: a factor equal to 3600 converting second level to hour level.
    parameter rou: the water density in kg/m^3.
    parameter lamb: the water vaporization in 2259.36 kJ/kg
    
    '''
    ratio = c/(rou*lamb)
    # Solar radiation approach
    if Approach=='Rs':
        ETd = (LE/Rsi)*Rsd*ratio
    # Evaporative fraction approach
    elif Approach=='EF':
        ETd = LE/(Rni-Gi)*(Rnd-Gd)*ratio
    # Ratio of net radiation to solar radiation approach
    elif Approach=='Rn-Rs':
        ETd = LE/(Rni-Gi)*Rni/Rsi*Rsd*ratio
    else:
        print('Please check the inputs.')
    
    return ETd