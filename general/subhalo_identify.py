#! /usr/bin/python
'''Scriptable module containing tools for identifying hosts and subhalos'''

from numpy import *

def subhalo_identify(halo_catalog):
    '''
    Identify host and subhalos and save this info to an IRATE halo
    catalog.

    This routine is addapted from a script written by Jose Onorbe. It
    identifies subhalos and adds a data set labeled "Hosts" in the given IRATE
    halo_catalog. Hosts = -1 for host halos that are not a subhalo of any other
    halo, all subhalos have Hosts != -1.

    Input:
     halo_catalog - A halo catalog from an IRATE file
    '''
    
    C = halo_catalog
    catalog_name = C.name
    print "\nStarting subhalo identification in catalog %s\n" % (catalog_name)
    
    # Easiest thing to do in this case is associate the x,y,z with the proper radius
    # I assume that the indices line up appropriately. I.e. for every center there is
    # a virial radius.
    
    centers = C["Center"][:]
    vradii  = C["Rvir"][:]
    if 'Mpc / h' in C['Center'].attrs['unitname']: vradii  = C["Rvir"][:]/1000.0
    vmass = C['Mvir'][:]

    # create some temporary host array
    tempHostArray = empty(len(centers),dtype=int)  
    #Now fill it with the -1 that we want:
    tempHostArray[:] = -1    


    #==========================
    # distance(center, center)
    #------------------------
    # 3D Distance equation; hard stuff.
    # Technically, I do not need the sqrt, and this could speed things up... depending on how the computers handle sqrt.
    def distance(a, b):
        return sqrt(pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2));


    #===============================
    # distanceNoSqrt(center, center)
    #-------------------------------
    # Calculates the square of a distance. The sqrt on some devices is an expensive function
    # and usually we're interested in just comparing magnitudes of distances and not being
    # incredibly precise. This is a possible optimization, which I don't really use in this
    # program, but I put it here for Shea to consider...
    def distanceNoSqrt(a,b):
        return pow(a[0]-b[0],2) + pow(a[1]-b[1],2) + pow(a[2]-b[2],2);

    #===============================
    # closer(center, center, center)
    #--------------------------------
    # calculates the distances of center j to center a, and center j to center b.
    # This tells us which center (a or b) is closer to j.
    # Notice that I only compare the squares of the distances. All that matters is the magnitudes
    # anyways.
    def closer(j, a, b):
        distja = distanceNoSqrt(centers[j], centers[a]);
        distjb = distanceNoSqrt(centers[j], centers[b]);

        if(distja < distjb):
            return a;
        elif (distjb < distja):
            return b;
        else:
            print("There is an equality for the distances, so this node qualifies for having two hosts...I'll return a - but this is weird. ");
            return a;

    # Inefficient, but easy, way of categorizing subhalos
    # It isn't quite enough to check that distances are less than some virial radius (though I do this), we need
    # to see whose radius is bigger (SHOULD IT BE BIGGER RVIR OR BIGGER MVIR?). So the condition should be that
    # if the distance is smaller than i's virial radius, and i's virial radius is larger than j's, then set i as j's host
    # in theory it's possible for a halo to seemingly have two hosts... what then? Well, we should pick the closer of the two as the host.
    # So, I append a check that says something like "If host is already assigned (i.e. not -1) then check the distance of this host
    # against the distance of the new host. Whichever one is closer gets to be the host.
    # Make sure that the halos don't mark themselves as their own host halo, especially since they will cross each other and their distance will be 0.
    print "Doing analysis...Now!";
    for i in range(len(centers)):
        myx,myy,myz = centers[i,0],centers[i,1],centers[i,2]  #xyz may be unnecessary, but it's more readable
        mymvir = vmass[i]
        distarray = sqrt(square(centers[:,0]-myx)+square(centers[:,1]-myy)+square(centers[:,2]-myz))    #So this is an array of all the distances, one of which (i) is identically zero
        #Now we want to pick halos using some numpy tricks, specifically where:
        subhalos = where(distarray<vradii[i])[0]  #For some reason, where gives the result as the [0]th entry in an array, hence the [0]
        #So subhalos is now a listing of all the halos that have centers within the virial radius of halo i

        #i.e. subhalos = all those halos selected by dist < vradii[i]

        #Unfortunately, it's not as straightforward as saying that all halos that have centers within the virial radius of halo i are subhalos of it
        #Instead, we need to check that all those halos have a virial mass smaller than the halo we're looking at
        for j in subhalos:
            if j == i:    continue   #We don't want to change the current halo

            #Let's check if it's already identified as a subhalo of something else
            if tempHostArray[j] != -1:  #Then it is a subhalo of something else
                #Let's check if what it's currently a subhalo of is more massive than the halo that we're currently looking at
                if vmass[tempHostArray[j]] > mymvir:
                    #Then the halo that we've assigned it to is more massive than the halo that we're anaylyzing
                    #Then that means (if subhalo defined only by mass) that the halo we're analyzing is probably a subhalo of the more massive one that it's already assigned to
                    #In that case, I'm going to choose to change it so that we have sub sub -> sub, which will eventually point to host
                    #We could check explicitely if the halo we're analyzing is a subhalo of the one we're removing it from, probably by first checking if it's host has been set then checking if it falls inside rvir of what we suspect is it's host, but I'm not sure if it's worth it.
                    #So, we're changing this one then:
                    tempHostArray[j] = i
                    continue    #Move on to the next subhalo
                else:
                    #Ok, so the halo that this subhalo is already assigned to is less massive than the halo that we're on
                    #So that means that the halo we're on is probably two generations above it, and it's already got it's true parent assigned
                    continue

            #now for the halos that haven't been identified as the subhalo of something else (this should be relatively straightforward)
            #Let's compare masses to determine the true host halo
            if mymvir > vmass[j]:   #So the halo that we're looking at is more massive
                tempHostArray[j] = i
            #else:   
                #then the halo that we've found to be within the virial radius is actually more massive than the halo we're looking at, so the halo we're looking at is a subhalo of halo j
                #but, it could be that j is two generations above it
                #With that in mind, and the knowledge that it'll get identified later anyways, let's just move on in this case.
                #Which means that we really don't even need this else statement

        #I think that actually does it, right?

        #Print statements are slow, so only do it ever 1000 or so
        if i%1000 == 0: #% is the modulus operator; i.e. the remainder
            print "%i has compared itself to everyone else."%i; # progress notifier

    print "Writing to HDF5 file...";
    if 'Hosts' in C.keys(): C['Hosts'][...] = tempHostArray
    else: hostsdset = C.create_dataset('Hosts',data=tempHostArray)
    print "Done.\n"
    
    return


## If used as a scritpt ##
if __name__ == "__main__":
    import sys
    import h5py
    if len(sys.argv) == 2:
        IrateFile = sys.argv[1]   
        snapnumber = -1
    elif len(sys.argv) == 3:
        IrateFile = sys.argv[1] 
        snapnumber =  sys.argv[2]
    else:   
        print "usage: python %s <IRATE_filename> [<snap_number>]" % __file__
        sys.exit(2)
          
    f = h5py.File(IrateFile, 'a')   
            
    if snapnumber == -1: snap_name = 'z=0'
    else: snap_name = snapnumber

    if IrateFile.find('snapshot')>0: # This is just to be compatible with old IRATE files
        C = f['Catalog']
    else:
        if snapnumber != -1:
            key = [key for key in f.keys() if snapnumber in key][0]
            snap = f[key]
        else:
            keys = [key for key in f.keys() if "Snapshot" in key]
            snap = f[sorted(keys)[-1]] # this will grab the Z=0 snapshot

        for key in snap.keys():
            if 'HaloCatalog' in key:
                catalog_name = key
                break
        C = snap[catalog_name]  
     
    subhalo_identify(C)    


