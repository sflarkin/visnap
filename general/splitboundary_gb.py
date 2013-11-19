#! /usr/bin/python
'''Scriptable module with tools to split the boundary particle group of gadget binary files'''

import sys, glob, os
import struct
from numpy import *
from visnap.general.read_gadget import read_gbin
#import pdb

def splitboundary_gb(infile, outfile, useStarGroup=False, AllBndryMasses=[], filedata=None):
    '''
    Read a gadget binary file with multple particle species in the boundary group
    and write a new gadget binary file placing the higher resolution particles of
    the boundary group into the disk, bulge and optionally the star groups.
    
     Input

      infile - name of gadget binary file to read

      outfile - name of the output splitted file to be written

      useStarGroup - If True the star group will be used when more
                     than 3 low-res particle species are found

      AllBndryMasses - If the file given is part of a multiple file gadget binary
                       the list AllBndryMasses tells the full range of masses in the
                       boundary group across all files.

      filedata - If the data of this file has already been read, you can provide
                 the filedata list of arrays containig the data to avoid calling
                 read_gbin() again.

    ''' 
    print "\nSplitting boundary group in file %s" % infile

    if filedata:
        [header,gasdata,halodata,diskdata,bulgedata,stardata,bndrydata,gasmasses] = filedata
    else:
        [header,gasdata,halodata,diskdata,bulgedata,stardata,bndrydata,gasmasses] = read_gbin(infile)
    

    if len(diskdata[0]) > 0 or len(bulgedata[0]) > 0:
        print "There are already particles in either disk or bulge.  Exiting."
        sys.exit()

    if len(bndrydata[0]) == 0:
        print "There are no boundary particles in this file, leaving as it is."
        os.system("cp "+infile+" "+outfile)
        return

    #My goal is to split the bndry group into disk, bulge, and optionally star based on masses

    thisfilebndrypos = bndrydata[0]
    thisfilebndryvel = bndrydata[1]
    thisfilebndryid = bndrydata[2]
    try:
        thisfilebndrymass = bndrydata[3]
    except IndexError:
        thisfilebndrymass = [header[1][5]]

    masses = []
    for i in range(len(thisfilebndrymass)):  #Loop through all the particles in the bndry group
        thismass = thisfilebndrymass[i]
        if thismass in masses:
            continue
        else:
            masses.append(thismass)
    masses.sort()    
    
    nmasses = len(masses)
    nAllBndryMasses = len(AllBndryMasses)

    if nAllBndryMasses == 0:
        AllBndryMasses = masses
        nAllBndryMasses = nmasses
    else:
        AllBndryMasses.sort()
        for m in masses:
            if m not in AllBndryMasses:
                print "The mass of one or more of the particle types in this file's boundary data in not in " \
                      " the AllBndryMasses list, this is odd! This file was not read when AllBndryMasses was set."
                sys.exit(1)
       
    mass_list = [0,0,0,0] # initialize [disk, bulge, star, bndry] masses to 0
    if useStarGroup:
        if nAllBndryMasses <= 4:
            mass_list[:nAllBndryMasses] = AllBndryMasses
        else:
            mass_list[:3] = AllBndryMasses[:3]
            mass_list[3] = AllBndryMasses[3:]
    else:    
        if nAllBndryMasses <= 2:
            mass_list[:nAllBndryMasses] = AllBndryMasses
        else:
            mass_list[:2] = AllBndryMasses[:2]
            mass_list[3] = AllBndryMasses[2:]
            
    diskmass, bulgemass, starmass, bndrymass = mass_list
    if not isscalar(bndrymass):
        if len(bndrymass) == 1:
            bndrymass = bndrymass[0]
        else:
            print "WARNING: The bndry group will have to be left with %d different low-res particle types" % (len(bndrymass))

    if starmass:
         print "WARNING: The star group will be used! If there is star formation in your simulation this won't work well."
         
    diskpos = thisfilebndrypos[thisfilebndrymass == diskmass]
    diskvel = thisfilebndryvel[thisfilebndrymass == diskmass]
    diskID = thisfilebndryid[thisfilebndrymass == diskmass]
    
    bulgepos = thisfilebndrypos[thisfilebndrymass == bulgemass]
    bulgevel = thisfilebndryvel[thisfilebndrymass == bulgemass]
    bulgeID = thisfilebndryid[thisfilebndrymass == bulgemass]
    
    starpos = thisfilebndrypos[thisfilebndrymass == starmass]
    starvel = thisfilebndryvel[thisfilebndrymass == starmass]
    starID = thisfilebndryid[thisfilebndrymass == starmass]

    if isscalar(bndrymass):
        bndrypos = thisfilebndrypos[thisfilebndrymass == bndrymass]
        bndryvel = thisfilebndryvel[thisfilebndrymass == bndrymass]
        bndryID = thisfilebndryid[thisfilebndrymass == bndrymass]
    else:
        bndrypos = empty([0,3],dtype=thisfilebndrypos.dtype)
        bndryvel = empty([0,3], dtype=thisfilebndryvel.dtype)
        bndryID = empty([0], dtype=thisfilebndryid.dtype )
        bndryMasses = empty([0],  dtype=thisfilebndrymass.dtype )
        for m in bndrymass:
            bndrypos = concatenate((bndrypos,thisfilebndrypos[thisfilebndrymass == m]))
            bndryvel = concatenate((bndryvel,thisfilebndryvel[thisfilebndrymass == m]))
            bndryID = concatenate((bndryID,thisfilebndryid[thisfilebndrymass == m]))
            bndryMasses = concatenate((bndryMasses,thisfilebndrymass[thisfilebndrymass == m]))
        bndrymass = 0    

    ndisk = len(diskpos)
    nbulge = len(bulgepos)
    nbndry = len(bndrypos)
    nstar = len(starpos)
    
    if nstar > 0 and header[0][4] > 0:
        print "Trying to place particles in the star group, but that group's not empty!  Please do this file by hand.  Exiting."
        sys.exit()
        
    elif nstar == 0:
        nstar = header[0][4]
        starmass = header[1][4]
        starpos = stardata[0]
        starvel = stardata[1]
        starID = stardata[2]

    if gasmasses:  print "There are {0} particles in the gas group with variable masses".format(header[0][0])
    else:        print "There are {0} particles in the gas group with a mass of {1}".format(header[0][0],header[1][0])
    print "There are {0} particles in the halo group with a mass of {1}".format(header[0][1],header[1][1])
    print "Placed {0} particles in the disk group with a mass of {1}".format(ndisk,diskmass)
    print "Placed {0} particles in the bulge group with a mass of {1}".format(nbulge,bulgemass)
    print "Placed {0} particles in the star group with a mass of {1}".format(nstar,starmass)

    if nbndry and bndrymass==0:
        print "Placed {0} particles in the bndry group with masses of {1}".format(nbndry,unique(bndryMasses))
    else:
        print "Placed {0} particles in the bndry group with a mass of {1}".format(nbndry,bndrymass)
                
    ntot = ndisk+nbulge+nbndry+header[0][1]+header[0][0]+nstar
    gas = header[0][0] > 0

    assert ndisk+nbulge+nbndry+nstar == header[0][5]      #Check that all the particles that were in the bndry group originally made it somewhere

    #pack the header and write it:
    
    f = open(outfile,'wb')
    f.write(struct.pack('<I',256))
    f.write(struct.pack('<6I',header[0][0],header[0][1],ndisk,nbulge,nstar,nbndry))    #Npart_file
    f.write(struct.pack('<6d',header[1][0],header[1][1],diskmass,bulgemass,starmass,bndrymass)) #Mass table 
    f.write(struct.pack('<d',header[2]))   #time
    f.write(struct.pack('<d',header[3]))   #z
    f.write(struct.pack('<i',header[4]))   #FlagSfr
    f.write(struct.pack('<i',header[5]))   #FlagFeedback
    f.write(struct.pack('<6i',header[6][0],header[6][1],ndisk,nbulge,nstar,nbndry)) #Npart_total
    f.write(struct.pack('<i',header[7]))   #FlagCooling
    f.write(struct.pack('<i',header[8]))   #NumFiles
    f.write(struct.pack('<d',header[9]))   #BoxSize
    f.write(struct.pack('<d',header[10]))  #Omega0
    f.write(struct.pack('<d',header[11]))  #OmegaL
    f.write(struct.pack('<d',header[12]))  #h
    f.write(struct.pack('<i',header[13]))  #FlagAge
    f.write(struct.pack('<i',header[14]))  #FlagMetals
    f.write(struct.pack('<6i',header[15][0],header[15][1],header[15][2],header[15][3],header[15][4],header[15][5]))    #NallHW -- might have to modify this for large sims
    f.write(struct.pack('<i',header[16]))  #flag_entr_ics

    header_bytes_left = 260 - f.tell()
    for i in range(header_bytes_left):
        f.write(struct.pack('<x'))    #Fill out the rest of the header with pad bytes
    assert f.tell() == 260
    f.write(struct.pack('<I',256))
    
    #Now pack the data using tostring, but I have to check that the dtypes of all my arrays are right
    gasdata[0] = gasdata[0].astype('f')
    gasdata[1] = gasdata[1].astype('f')
    gasdata[2] = gasdata[2].astype('I')
    if gas:
        gasdata[3] = gasdata[3].astype('f')      #If gasmasses, this is the mass block, otherwise, it's the internal energy
        if gasmasses:
            gasdata[4] = gasdata[4].astype('f')  #This is the internal energy if it exists, but it only exists if gasmasses is true
    
    halodata[0] = halodata[0].astype('f')
    halodata[1] = halodata[1].astype('f')
    halodata[2] = halodata[2].astype('I')
    try:
        halodata[3] = halodata[3].astype('f')
    except IndexError:
        print "There is no mass block for the halo group, as expected."

    bndrypos = bndrypos.astype('f')
    bndryvel = bndryvel.astype('f')
    bndryID = bndryID.astype('I')
    
    bulgepos = bulgepos.astype('f')
    bulgevel = bulgevel.astype('f')
    bulgeID = bulgeID.astype('I')
    
    diskpos = diskpos.astype('f')
    diskvel = diskvel.astype('f')
    diskID = diskID.astype('I')
    
    starpos = starpos.astype('f')
    starvel = starvel.astype('f')
    starID = starID.astype('I')
    try:
        stardata[3] = stardata[3].astype('f')
    except IndexError:
        print "There is no mass block for star particles, as expected."
        
    #Now I have to write it all to the file
    print "Writing coordinates."
    f.write(struct.pack('<I',12*ntot))
    f.write(gasdata[0].tostring())
    f.write(halodata[0].tostring())
    f.write(diskpos.tostring())
    f.write(bulgepos.tostring())
    f.write(starpos.tostring())
    f.write(bndrypos.tostring())
    f.write(struct.pack('<I',12*ntot))
    
    
    print "Writing velocities."
    f.write(struct.pack('<I',12*ntot))
    f.write(gasdata[1].tostring())
    f.write(halodata[1].tostring())
    f.write(diskvel.tostring())
    f.write(bulgevel.tostring())
    f.write(starvel.tostring())
    f.write(bndryvel.tostring())
    f.write(struct.pack('<I',12*ntot))
    
    
    print "Writing particles IDs."
    f.write(struct.pack('<I',4*ntot))
    f.write(gasdata[2].tostring())
    f.write(halodata[2].tostring())
    f.write(diskID.tostring())
    f.write(bulgeID.tostring())
    f.write(starID.tostring())
    f.write(bndryID.tostring())
    f.write(struct.pack('<I',4*ntot))
    
    nmass = 0
    if len(stardata) > 3:
        nmass = nmass + len(stardata[3])
    if len(halodata) > 3:
        nmass = nmass + len(halodata[3])
    if gasmasses:
        nmass = nmass + len(gasdata[3])
    if nbndry and bndrymass==0:
        nmass = nmass + len(bndryMasses)
        
    if gasmasses or (header[1][1] == 0 and header[0][1] > 0) or (header[1][4] == 0 and header[0][4] > 0) or (nbndry and bndrymass==0):
        #Then I have to do masses too
        f.write(struct.pack('<I',4*nmass))
        if gasmasses:
            print "Writing variable masses for gas particles."
            f.write(gasdata[3].tostring())
        if len(halodata) > 3:
            print "Writing variable masses for halo particles."
            f.write(halodata[3].tostring())
        if len(stardata) > 3:
            print "Writing variable masses for star particles."
            f.write(stardata[3].tostring())
        if nbndry and bndrymass==0:
            print "Writing variable masses for boundary particles."
            f.write(bndryMasses.tostring())
        f.write(struct.pack('<I',4*nmass))
    else:
        print "There are no particles with variable masses being written."
    
    
    if gas:
        print "Writing gas internal energy."
        #Then I have to do gas specific stuff
        f.write(struct.pack('<I',4*len(gasdata[0])))
        if gasmasses:   f.write(gasdata[4].tostring())
        else:         f.write(gasdata[3].tostring())
        f.write(struct.pack('<I',4*len(gasdata[0])))
    
    f.close()
    print 'wrote '+outfile


def splitboundary_gbs(inbase, outbase, useStarGroup=False):
    '''
    Read a multiple file gadget binary with multple particle species in the
    boundary group and write a new multiple file gadget binary placing the
    higher resolution particles of the boundary group into the disk, bulge
    and optionally the star groups.
    
     Input

      inbase - basename of gadget binary files to read

      outbase - basename of the output splitted files to be written

      useStarGroup - If True the star group will be used when more
                     than 3 low-res particle species are found
    '''
   
    infiles = glob.glob(inbase+'*')
    skipFiles = []
    AllBndryMasses=[]
    FilesData = []

    print "Reading all files to find the range of masses in the boundary group"
    for infile in infiles:
        [header,gasdata,halodata,diskdata,bulgedata,stardata,bndrydata,gasmasses] = thisfiledata = read_gbin(infile)
        FilesData.append(thisfiledata)
        
        if len(diskdata[0]) > 0 or len(bulgedata[0]) > 0:
            print "There are already particles in either the disk or bulge of file %s. Exiting." % infile
            sys.exit()

        if len(bndrydata[0]) == 0:
            ext = infile.split('.')[-1]
            outfile = outbase+'.'+ext
            print "There are no boundary particles in file %s, copying as it is to file %s" % (infile, outfile)
            os.system("cp "+infile+" "+outfile)
            skipFiles.append(infile)
            continue

        try:
            thisfilebndrymass = bndrydata[3]
        except IndexError:
            thisfilebndrymass = [header[1][5]]
            
        for i in range(len(thisfilebndrymass)):  #Loop through all the particles in the bndry group
            thismass = thisfilebndrymass[i]
            if thismass in AllBndryMasses:
                continue
            else:
                AllBndryMasses.append(thismass)
            
    AllBndryMasses.sort()        

    print "\nSplitting files"        
    for i,infile in enumerate(infiles):
        if infile in skipFiles:
            continue
        else:
            ext = infile.split('.')[-1]
            outfile = outbase+'.'+ext
            splitboundary_gb(infile, outfile, useStarGroup=useStarGroup,
                     AllBndryMasses=AllBndryMasses, filedata=FilesData[i])
    


if __name__ == "__main__":

    if len(sys.argv) < 3 or len(sys.argv) > 4 :
        print "Usage: python {0} <input file(base)> <output file(base)> [options]".format(sys.argv[0])
        print "Options:"
        print " -s  Use star group"
        sys.exit()
              
    infiles = glob.glob(sys.argv[1]+'*')
    
    if '-s' in sys.argv:
        useStarGroup=True
    else:
        useStarGroup=False 
    
    if len(infiles) == 1:
        splitboundary_gb(infiles[0], sys.argv[2], useStarGroup=useStarGroup)
    else:
        splitboundary_gbs(sys.argv[1], sys.argv[2], useStarGroup=useStarGroup)


