'''Module containing tools that translate simulation file names to simulation properties'''

def hyades(filename):
    '''
    Exctract simulation parameters from hyades simulations filenames

    Input:
     filename - A simulation filename following the hyades simulations name
     convention

    Output: 
     A dictionary containing the following simulation properties.
      dmName - CDM, SIDM, WDM or DDM
      Lbox - Box size
      Nbox - Number of particles 
      mpdm - Mass of dark matter particles
      epsilon - Force resolution
      zoom_id - Halo id for Zooms 
      sigma_m - DM self-interaction cross section [cm^2/g]
      hsi - SIDM smoothing factor 
      wdmMass - Warm DM particle mass
      DDMtau -  Decaying DM time scale
      DDMvk - Decaying DM kick velocity 
    '''
    name = filename.split('/')[-1]
    nameSplit = name.split('_')
    if nameSplit[0] == 'GID':  dmName = 'SIDM'
    if nameSplit[0] == 'GVD':  dmName = 'CDM'
    if nameSplit[0] == 'GWD':  dmName = 'WDM'
    if nameSplit[0] == 'GDD':  dmName = 'DDM'
    if nameSplit[1].isdigit():
        if nameSplit[1][0]==0: Lbox = nameSplit[1][1:]
        else:  Lbox = nameSplit[1]
        Nbox = 2**float(nameSplit[2])
        mpdm = 6.88*(10**7)*(512.0/50.0)**3*(float(Lbox)/Nbox)**3
        Nbox = str(Nbox)
        if  int(nameSplit[3][0]) == 1:  epsilon = float(nameSplit[3][1:])
        elif int(nameSplit[3][0]) == 0:  epsilon = float(nameSplit[3][1:])/1000.0
        else: 
            print "The softening identifier in the file name %s does't start with 0 nor 1, I don't know what that means." % ifile
            sys.exit(1)
        zoom_id = nameSplit[4].lstrip('0')
        if dmName == 'WDM':
            wdmMass = nameSplit[5]
        else:
            wdmMass = None
        if dmName == 'DDM':
            DDMtau, DDMvk = float(nameSplit[5][1:]), float(nameSplit[6][1:])
        else: 
            DDMtau, DDMvk = None,None
        if dmName == 'SIDM':
            sigma_m =  float(nameSplit[5][1] + '.' + nameSplit[5][2:])
            if nameSplit[6][0] == 'h': hsi = float(nameSplit[6][1:])
            else: hsi = None
        else:
            sigma_m = None
            hsi = None
    else:    
        if nameSplit[1] == 'D08':  Lbox, Nbox, mpdm = '25', '256',6.88e7
        if nameSplit[1] == 'D09':  Lbox, Nbox, mpdm = '25', '512',8.59e6
        if nameSplit[1] == 'F08':  Lbox, Nbox, mpdm = '50', '256',5.50e8
        if nameSplit[1] == 'F09':  Lbox, Nbox, mpdm = '50', '512',6.88e7
        if nameSplit[1] == 'Z11':  Lbox, Nbox, mpdm = '50', '2048',1.07e6
        if nameSplit[1] == 'Z12':  Lbox, Nbox, mpdm = '50', '4096',1.34e5 
        if dmName == 'WDM': 
            wdmMass = nameSplit[4]
        else:
            wdmMass = None
        if dmName == 'DDM':
            DDMtau, DDMvk = float(nameSplit[4][1:]), float(nameSplit[5][1:])
        else: 
            DDMtau, DDMvk = None,None
        if dmName == 'SIDM':    
            sigma_m =  float(nameSplit[4][1] + '.' + nameSplit[4][2:])
            if nameSplit[5][0] == 'h': hsi = float(nameSplit[5][1:])
            else: hsi = None
        else:    
            sigma_m = None
            hsi = None
        if  int(nameSplit[2][0]) == 1:  epsilon = float(nameSplit[2][1:])
        elif int(nameSplit[2][0]) == 0:  epsilon = float(nameSplit[2][1:])/1000.0
        else: 
            print "The softening identifier in the file name %s does't start with 0 nor 1, I don't know what that means." % ifile
            sys.exit(1)
        zoom_id = nameSplit[3].lstrip('0')  
    mpdmS = ('%.1e' % mpdm).split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')   
    nameDict = {'dmName': dmName, 'Lbox': Lbox, 'Nbox': Nbox, 'mpdm': mpdm, 
                'epsilon': epsilon, 'zoom_id': zoom_id, 'sigma_m': sigma_m, 'hsi': hsi, 
                'wdmMass': wdmMass, 'DDMtau': DDMtau, 'DDMvk': DDMvk}
    return nameDict

