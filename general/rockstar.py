'''Module containing tools useful to handle Rockstar catalogs and files'''

from numpy import *

def read_rockstar_header(rockstar_file):
    '''
    Read rockstar ASCII headers
    
    Input:
     rockstar_file - An ASCII rockstar file, either halo catalog or merger tree
     
    Output:
      comments, cosmo, unitnames, unitcgs - Dictionaries ready to use to write
                                            to an IRATE file as attributes 
    '''
    f = open(rockstar_file)
    linen = 0
    comments = []
    units = []
    cosmo = []
    try:
        for line in f:
            if  (line[0] == "#") and (':' not in line) and  (';' not in line) and  ('=' not in line): 
                if linen == 0:
                    columns = line.lstrip('#')
                    columns = columns.split()
                    linen = 1
                continue
            elif line[0] == "#":
                thisl = line.lstrip('#').rstrip('\n')
                if not thisl.startswith('Notes'):
                    thisl = thisl.split(';')
                else: thisl = [thisl]    
                for j in range(len(thisl)):
                    thiscom = thisl[j]
                    thiscom = thiscom.strip()
                    if not thiscom.startswith('Units'):
                        #So, not a units line
                        if thiscom.startswith('Om =') or thiscom.startswith('Ol =') or thiscom.startswith('h ='):
                            #Then it's cosmology, and I want to save that separately
                            thiscom = thiscom.split('=')
                            thiscom[0] = thiscom[0].strip()
                            thiscom[-1] = float(thiscom[-1])
                            if thiscom not in cosmo:        #skip it if it's already in there, which would come up if files from multiple processors are catted together
                                cosmo.append(thiscom)
                        else:
                            thiscom = thiscom.split(':')
                            if len(thiscom) == 1:   #Then it didn't split on a :. Either it has to be split on =, or not a comment 
                                if '=' in thiscom[0]:
                                    thiscom = thiscom[0]
                                    thiscom = thiscom.split('=')
                                    for k in range(len(thiscom)):
                                        thiscom[k] = thiscom[k].strip()        
                            if len(thiscom)==2:            
                                if (thiscom not in comments) and (thiscom[1] != ''):  #again, skip if it's already been read, or if empty
                                    comments.append(thiscom)
                            else:
                                print 'WARNING: The comment %s was not supported. It will not be included in the attributes'%(thiscom)        
                    else:
                        #What to do if it's a line about units
                        thiscom = thiscom.lstrip('Units: ').split(' in ')
                        if len(thiscom) == 2:
                            mycom = [thiscom[0],thiscom[1]]
                        #elif len(thiscom) < 2:
                        #    thiscom.split(' are ') # How to treat dimensionless units lines? for now we simply dont include them 
                        #    mycom = [thiscom[0],thiscom[1]]
                        if mycom not in units:          #skip if it's already been read in
                            units.append(mycom)
                        
            else:  break
    finally: f.close()

    for i in range(len(units)):
        #Figure out the units for each thing
        if  "Masses" in units[i][0]:
            munitname = units[i][1]
        elif "Positions" in  units[i][0]:
            posunitname = units[i][1]
        elif "Velocities" in units[i][0]:
            velunitname = units[i][1]
        elif "Radii" in units[i][0]:
            radunitname = units[i][1]
        elif "Angular Momenta" in units[i][0]:
            angunitname = units[i][1]
        elif "Total energy" in units[i][0]:
            enunitname = units[i][1]
        else:
            print "There is a unit in the file that wasn't planned for.  Note that this may mean that certain datasets won't have units, which will make this an invalid IRATE file."

    #Now let's make the arrays for the default values and just spit out a warning if it's not the default
    if munitname == "Msun / h":  munitcgs = array([1.98892e33,-1,0])
    else:   
        print "Unrecognized unit for mass--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        munitcgs = None
    
    if posunitname == "Mpc / h (comoving)": posunitcgs = array([3.08568025e24,-1,1])
    else:   
        print "Unrecognized unit for position--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        posunitcgs = None
    
    if velunitname == "km / s (physical, peculiar)" or velunitname == "km / s (physical)":  velunitcgs = array([1e5,0,0])
    else:
        print "Unrecognized unit for velocity--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        velunitcgs = None
    
    if radunitname == "kpc / h (comoving)":  radunitcgs = array([3.08568025e21,-1,1])
    else:
        print "Unrecognized unit for radii--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        radunitcgs = None
        
    if angunitname == "(Msun/h) * (Mpc/h) * km/s (physical)":  angunitcgs = array([1.98892e33*3.08568025e24*1e5,-2,0])
    else:
        print "Unrecognized unit for angular momenta--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        angunitcgs = None
    
    if enunitname == "(Msun/h)*(km/s)^2 (physical)":  enunitcgs = array([1.98892e33*1e5*1e5,-1,0])
    else:
        print "Unrecognized unit for total energy--unitname will be saved, but no numerical conversion factors will be supplied, so this won't be a valid IRATE file."
        enunitcgs = None

    unitnames = {'Mass': munitname, 'Position': posunitname, 'Velocity': velunitname, 'Radius': radunitname, 'AngularM': angunitname, 'Energy': enunitname }
    unitcgs =  {'Mass': munitcgs, 'Position': posunitcgs, 'Velocity': velunitcgs, 'Radius': radunitcgs, 'AngularM': angunitcgs, 'Energy': enunitcgs }

    return comments, cosmo, unitnames, unitcgs


