#!/usr/bin/env python

import h5py,os
from optparse import OptionParser
from numpy import *
import sys

usage = "usage: %prog [options] <input-file> <snapshot-number> <output-file> <xmin> <ymin> <zmin> <dx> <dy> <dz> \n or %prog [options] <input-file> <snapshot-number> <output-file>  <xcenter> <ycenter> <zcenter> <R> "
parser =  OptionParser(usage=usage)

parser.add_option("-n","--data-name",dest="dname",help="The name to be given to the new group containing the data for particles within the given region. Default = %default",default='Dark_Halo_Region')


(ops,args) = parser.parse_args()
dname = ops.dname

try:
    args
except NameError:
    args = sys.argv

if '--' in args: args.remove('--')
sphere = 0

if len(args) == 9:
    inname,snapnum,outname = args[0:3]
    xmin,ymin,zmin,dx,dy,dz = asfarray(args[3:])
elif len(args) == 7:
    inname,snapnum,outname = args[0:3]
    xcenter, ycenter, zcenter, R =  asfarray(args[3:])
    sphere =1
else:
    parser.print_help()
    sys.exit(1337)

try:
    snapnum = int(snapnum)
    snapnumS = str(snapnum).zfill(5)
except ValueError:
    parser.print_help()
    print "The value you provided for the snapshot number is not an integer.  Please try again."
    sys.exit(1337)

#open input irate file and read headers

inirate = h5py.File(inname,'r')
insimprops = inirate['SimulationProperties']
incosmo = inirate['Cosmology']

#Get all the data
try:
    insnap = inirate['Snapshot'+snapnumS]
    inpdata = insnap['ParticleData']
    inhalo = inpdata['Dark_Halo']
except KeyError:
    print "There is no group /Snapshot/%s/ParticleData/Dark_Halo.  Please check your file format and the snapshot number you provided and try again." % snapnumS
    sys.exit(1337)

dpos = inhalo['Position'][...]
dvel = inhalo['Velocity'][...]
dID = inhalo['ID'][...]
dmass = inhalo['Mass'][...]

print 'x_max,xmin = ', dpos[:,0].max(), dpos[:,0].min()
print 'y_max,ymin = ', dpos[:,1].max(), dpos[:,1].min()
print 'z_max,z_min = ', dpos[:,2].max(), dpos[:,2].min()

#find the particles that you want to slect for rewriting into the new file

if sphere==1:
    print "slecting spherical region centered at [%g,%g,%g] with Radius = %g" % (xcenter,ycenter,zcenter,R)
    iwant = argwhere(( abs(dpos[:,0]-xcenter) < R ) & ( abs(dpos[:,1]-ycenter) < R) & (abs(dpos[:,2]-zcenter) < R))[:,0]

if sphere == 0:
    print "slecting rectangular region with [xmin,ymin,zmin] = [%g,%g,%g] and [dx,dy,dz] = [%g,%g,%g]" % (xmin,ymin,zmin,dx,dy,dz)
    iwant = argwhere( (dpos[:,0] > xmin) & (dpos[:,0] < xmin + dx) & (dpos[:,1] > ymin) & (dpos[:,1] < ymin + dy) &  (dpos[:,2] > zmin) & (dpos[:,2] < zmin + dz) )[:,0]

print "\nSelected %d particles within the given region\n"%iwant.size

#Now we have rewrite this selection into an irate file

#start by rewriting the things that don't change i.e. cosmology and simprops

if os.path.isfile(outname):
    outirate = h5py.File(outname,'a')
    print "Opening {0} to add data to it.".format(outname)
    if "Cosmology" not in outirate.keys():
        #add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=True)
        outcosmo = outirate.create_group('Cosmology')
        for key in incosmo.attrs.keys(): outcosmo.attrs[key] = incosmo.attrs[key]
    #else:
    #    try:
    #        add_cosmology(outname,omegaM=omegaM,omegaL=omegaL,h=hubble,s8=s8,ns=ns,omegaB=omegaB,update=False)
    #    except KeyError:    #Then ns, s8, or omegaB was passed, but wasn't in the file before (in which case I want to add it)
    #        add_cosmology(outname,s8=s8,ns=ns,omegaB=omegaB,update=True)
    if "SimulationProperties" not in outirate.keys():
        outsimprops = outirate.create_group( "SimulationProperties")
        for key in insimprops.attrs.keys(): outsimprops.attrs[key] = insimprops.attrs[key]
else:
    print "Creating new IRATE file {0}".format(outname)
    outirate = h5py.File(outname,'w')
    outcosmo = outirate.create_group('Cosmology')
    outsimprops = outirate.create_group( "SimulationProperties") 
    for key in incosmo.attrs.keys(): outcosmo.attrs[key] = incosmo.attrs[key] 
    for key in insimprops.attrs.keys(): outsimprops.attrs[key] = insimprops.attrs[key] 

#Now lets write our particle selection

try:
    outsnap = outirate.create_group("Snapshot"+snapnumS)
    for key in insnap.attrs.keys(): outsnap.attrs[key] = insnap.attrs[key]
    print "Created new group for Snapshot "+snapnumS
except ValueError:
    outsnap = outirate['Snapshot'+snapnumS]
    print "Adding data to existing group Snapshot"+snapnumS

try:
    outpdata = outsnap.create_group('ParticleData')
    #Add units to pdata
    for key in inpdata.attrs.keys(): outpdata.attrs[key] = inpdata.attrs[key]
except ValueError:
    outpdata = outsnap['ParticleData']
    #If 'ParticleData' already exists, check that the units that I'm adding match up (if it has units)
    try:
        #Check that the units that are there match what I have
        msg = "\n"
        exitflag = False
        if outpdata.attrs['Positionunitname'] != inpdata.attrs['Positionunitname'] or (outpdata.attrs['Positionunitcgs'] != inpdata.attrs['Positionunitcgs']).all():
            msg = msg + "Length units in existing file are different than those provided.\n"
            exitflag = True
        if outpdata.attrs['Velocityunitname'] != inpdata.attrs['Velocityunitname']  or (outpdata.attrs['Velocityunitcgs'] != inpdata.attrs['Velocityunitcgs']).all():
            msg = msg + "Velocity units in existing file are different than those provided.\n"
            exitflag = True
        if outpdata.attrs['Massunitname'] !=  inpdata.attrs['Massunitname'] or (outpdata.attrs['Massunitcgs'] != inpdata.attrs['Massunitcgs'] ).all():
            msg = msg + "Mass units in existing file are different than those provided.\n"
            exitflag = True
        if exitflag:
            raise ValueError(msg)
    except KeyError:
    #add the units, since they're not there 
        outpdata.attrs['Positionunitname'] = inpdata.attrs['Positionunitname']
        outpdata.attrs['Positionunitcgs'] = inpdata.attrs['Positionunitcgs'] 
        outpdata.attrs['Velocityunitname'] = inpdata.attrs['Velocityunitname']
        outpdata.attrs['Velocityunitcgs'] = inpdata.attrs['Velocityunitcgs']
        outpdata.attrs['Massunitname'] = inpdata.attrs['Massunitname']
        outpdata.attrs['Massunitcgs'] =  inpdata.attrs['Massunitcgs']

exitflag = False
capbet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

if dname[0] not in capbet:  dname = dname.capitalize()
try:
    outhalo = outpdata.create_group(dname)
    outhalo.attrs['Npart'] = iwant.size
    print "Saving halo data under /Snapshot{0:05}/ParticleData/{1}".format(snapnum,dname)
except ValueError:
    msg = msg + "The group /Snapshot{0:05}/ParticleData/{1} already exists. Please specify a different location to save halo data.\n".format(snapnum,dname)
    exitflag = True

if exitflag:
    raise ValueError(msg)   

outhalo.create_dataset('Position',data=dpos[iwant])
outhalo.create_dataset('Velocity',data=dvel[iwant])
outhalo.create_dataset('ID',data=dID[iwant])
outhalo.create_dataset('Mass',data=dmass[iwant])

inirate.close()
outirate.close()

