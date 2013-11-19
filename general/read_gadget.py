'''Module containing tools for reading gadget files'''

import struct
from numpy import fromstring

def read_gbin(fname):
    '''
    Read gadget binary file

    Input:

     fname - name of the gadget binary file to read

    Output:
    
     [header,gas,halo,disk,bulge,star,bndry,gasmass] - arrays containing the gadget
                                                       data structures

    '''
    
    print "\nOpening "+fname    
    f = open(fname,'rb')

    #First read the header so I know how many particles of each type I have
    header_size = struct.unpack('<I',f.read(4))[0]

    #number of particles of each type in this file
    nfile = struct.unpack('<6I',f.read(24)) #Number of particles in this file

    masstable = struct.unpack('<6d',f.read(48))  #masses of the particle groups
        
    a = struct.unpack('<d',f.read(8))[0]        #expansion factor
    z = struct.unpack('<d',f.read(8))[0]        #redshift

    flag_sfr = struct.unpack('<i',f.read(4))[0] #star formation included?
    flag_feed = struct.unpack('<i',f.read(4))[0] #feedback included?

    ntot = struct.unpack('<6i',f.read(24))      #total number of particles in the simulation (= nfile if numfiles == 1)
        
    flag_cool = struct.unpack('<i',f.read(4))[0]  #cooling included?
    numfiles = struct.unpack('<i',f.read(4))[0]   #number of files in each snapshot
    boxsize = struct.unpack('<d',f.read(8))[0] #Size of the box, if periodic
    omega0 = struct.unpack('<d',f.read(8))[0]  #matter density at z = 0
    omegaL = struct.unpack('<d',f.read(8))[0]  #vacuum energy density at z = 0
    h = struct.unpack('<d',f.read(8))[0] #hubble parameter in units of 100 km/s/Mpc
    flag_age = struct.unpack('<i',f.read(4))[0]  #stellar age included?
    flag_metals = struct.unpack('<i',f.read(4))[0]  #use metals?
    nhighword = struct.unpack('<6i',f.read(24))   #contains the most significant word of 64-bit particle numbers (if npart > 2^32)

    flag_entropy = struct.unpack('<i',f.read(4))[0] #entropy instead of thermal energy in initial conditions?

    f.seek(264,0)   #Moves to the end of the header (and block that tells you size of header)

    header = [nfile,masstable,a,z,flag_sfr,flag_feed,ntot,flag_cool,numfiles,boxsize,omega0,omegaL,h,flag_age,flag_metals,nhighword,flag_entropy]
    
    gas,halo,disk,bulge,star,bndry = [],[],[],[],[],[]
    ngas = nfile[0]
    nhalo = nfile[1]
    ndisk = nfile[2]
    nbulge = nfile[3]
    nstar = nfile[4]
    nbndry = nfile[5]
    
    print "Reading coordinates"
    coord_size = struct.unpack('<I',f.read(4))[0]    
    gas.append(fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    halo.append(fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    disk.append(fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    bulge.append(fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    star.append(fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    bndry.append(fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))
    print "Read coordinates for {0} of {1} particles".format(ngas+nhalo+ndisk+nbulge+nstar+nbndry,sum(nfile))
    if struct.unpack('<I',f.read(4))[0] != coord_size:
        raise StandardError("The block size at the end of the coordinate block doesn't match that at the beginning.  This is an issue.")

    #next up is velocities.  pretty identical to the coordinates.
    vel_size = struct.unpack('<I',f.read(4))[0]
    print "Reading velocities"
    gas.append(fromstring(f.read(12*ngas),dtype='f').reshape((-1,3)))
    halo.append(fromstring(f.read(12*nhalo),dtype='f').reshape((-1,3)))
    disk.append(fromstring(f.read(12*ndisk),dtype='f').reshape((-1,3)))
    bulge.append(fromstring(f.read(12*nbulge),dtype='f').reshape((-1,3)))
    star.append(fromstring(f.read(12*nstar),dtype='f').reshape((-1,3)))
    bndry.append(fromstring(f.read(12*nbndry),dtype='f').reshape((-1,3)))
    print "Read velocities for {0} of {1} particles".format(ngas+nhalo+ndisk+nbulge+nstar+nbndry,sum(nfile))
    if struct.unpack('<I',f.read(4))[0] != vel_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the velocity block doesn't match that at the beginning.  This is an issue.")

    #Next up are the particle IDs, which are unsigned integers.
    id_size = struct.unpack('<I',f.read(4))[0]
    print "Reading particle IDs"
    gas.append(fromstring(f.read(4*ngas),dtype='I'))
    halo.append(fromstring(f.read(4*nhalo),dtype='I'))
    disk.append(fromstring(f.read(4*ndisk),dtype='I'))
    bulge.append(fromstring(f.read(4*nbulge),dtype='I'))
    star.append(fromstring(f.read(4*nstar),dtype='I'))
    bndry.append(fromstring(f.read(4*nbndry),dtype='I'))
    print "Read IDs for {0} of {1} particles".format(ngas+nhalo+ndisk+nbulge+nstar+nbndry,sum(nfile))
    if struct.unpack('<I',f.read(4))[0] != id_size:   #And read the size of the block again.
        raise StandardError("The block size at the end of the IDs block doesn't match that at the beginning.  This is an issue.")

    #Now I have to do the mass block.  Do this by checking if there are particles in a block that have a zero in the mass table.    
    if ((ngas > 0 and masstable[0] == 0) or (nhalo > 0 and masstable[1] == 0)
        or (ndisk > 0 and masstable[2] == 0) or (nbulge > 0 and masstable[3] == 0)
        or (nstar > 0 and masstable[4] == 0) or (nbndry > 0 and masstable[5] == 0)):
        #In other words, only read the size of the mass block if there is a mass block (for any of the groups)
        mass_size = struct.unpack('<I',f.read(4))[0]
        
        if ngas > 0 and masstable[0] == 0:    #There are particles in the group, but their masses aren't in the header (so they must be in the file)
            print "Reading variable masses for gas group"
            gas.append(fromstring(f.read(4*ngas),dtype='f'))
            gasmasses = True
        else:
            gasmasses = False
            
        if nhalo > 0 and masstable[1] == 0:
            print "Reading variable masses for halo group"
            halo.append(fromstring(f.read(4*nhalo),dtype='f'))
 
        if ndisk > 0 and masstable[2] == 0:
            print "Reading variable masses for disk group"
            disk.append(fromstring(f.read(4*ndisk),dtype='f'))
 
        if nbulge > 0 and masstable[3] == 0:
            print "Reading variable masses for bulge group"
            bulge.append(fromstring(f.read(4*nbulge),dtype='f'))
 
        if nstar > 0 and masstable[4] == 0:
            print "Reading variable masses for star group"
            star.append(fromstring(f.read(4*nstar),dtype='f'))
 
        if nbndry > 0 and masstable[5] == 0:
            print "Reading variable masses for boundary group"
            bndry.append(fromstring(f.read(4*nbndry),dtype='f'))
 
        if struct.unpack('<I',f.read(4))[0] != mass_size:   #And read the size of the block again.
            raise StandardError("The block size at the end of the mass block doesn't match that at the beginning.  This is an issue.")
    else:
        gasmasses = False
        

    if ngas > 0:
        print "Reading gas specific data."
        
        #Put all this inside the if statement because I don't want it reading block sizes for gas blocks if there is no gas data (cause then there will be no block size and I will accidentally read into other data or past the end of the file.)
        
        #Internal energy:
        u_size = struct.unpack('<I',f.read(4))[0]
        gas.append(fromstring(f.read(4*ngas),dtype='f'))
        if struct.unpack('<I',f.read(4))[0] != u_size:
            raise StandardError("The block size at the end of the internal energy block doesn't match that at the beginning.  This is an issue.")
        print "Read gas internal energy; not attempting to read gas density or smoothing length because this is assumed to be an IC file"
        
        
    current_pos = f.tell()
    f.seek(0,2) #Jump to the end of the file
    if current_pos == f.tell():
        print "Read "+fname
    else:
        print "Completed reading "+fname+" but there remain {0} bytes at the end of the file unread.".format(f.tell()-current_pos)
    f.close()

    return [header,gas,halo,disk,bulge,star,bndry,gasmasses]


