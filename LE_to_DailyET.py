def LE_to_DailyET(LE,Approach='Rs',
                  Rsi=500,Rsd=500,
                  Rni=500,Rnd=500,
                  Gi=500,Gd=500,
                  c=1000,rou=997,lamb=2259.36):
    '''
    parameter LE: the instantent latent heat flux in W/m2.
    parameter Approach: this one is string and it is restricted in 'Rs', 'EF', and 'Rn-Rs'.
    
    parameter Rsi: the instantent solar radiation in W/m2.
    parameter Rsd: the summation of the instantent solar radiation in W/m2.
    parameter Rni: the instantent net radiation in W/m2.
    parameter Rnd: the summation of the instantent net radiation in W/m2.
    parameter Gi: the instantent soil heat flux in W/m2.
    parameter Gd: the summation of the instantent soil heat flux in W/m2.
    
    parameter c: a factor equal to 1000 converting meter to millimeter.
    parameter rou: the water density in kg/m3.
    parameter lamb: the water vaporization in MJ/kg. The heat of vaporization of water is around 540 cal/g(https://www.google.com/search?q=the+heat+of+vaporization+of+water&rlz=1C1CHBF_enUS932US932&oq=the+heat+of+vaporization+of+water&aqs=chrome..69i57j0l6j0i395l3.6550j1j4&sourceid=chrome&ie=UTF-8), which is around 2259.36MJ/kg.
    
    '''
    ratio = c/(rou*lamb)
    # Solar radiation approach
    if Approach=='Rs':
        ETd = LE/Rsi*Rsd*ratio
    # Evaporative fraction approach
    elif Approach=='EF':
        ETd = LE/(Rni-Gi)*(Rnd-Gd)*ratio
    # Ratio of net radiation to solar radiation approach
    elif Approach=='Rn-Rs':
        ETd = LE/(Rni-Gi)*Rni/Rsi*Rsd*ratio
    else:
        print('Please check the inputs.')
    
    return ETd