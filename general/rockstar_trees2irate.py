#! /usr/bin/python
'''Scriptable module with tools to write rockstar trees into an IRATE file'''
import pdb

def read_trees(fname):
    '''
    Read rockstar merger trees and write them to an IRATE file 
    
    Input:
     fname - File with rockstar trees

    Output:
     colheads,trees,treenums,comments - As required for save_to_hdf5()
    
    '''
    import os,sys
    try:
        from progressbar import ProgressBar,Bar,Percentage,ETA
        dopbar = True
    except ImportError:
        dopbar = False    
#    if os.path.isfile('tmp.txt'):
#        print "There's a file in the current directory named 'tmp.txt', which would get overwritten by this script.  Please rename or remove that file."
#        sys.exit(1337)
    from numpy import loadtxt
    from StringIO import StringIO
    print "Opening {0}".format(fname)
    f = open(fname,'r')
    thistree = ""
    comments = []
    treenums = []
    trees = []
    linen = 0
    ex = 0
    startedtrees = False
    for line in f:
        if linen == 0:
            linen = linen + 1
            colheads = line.strip('#')
            colheads = colheads.split()
            for i in range(len(colheads)):
                colheads[i] = colheads[i].rstrip('(0123456789)')
        elif line.startswith('#tree'):
            startedtrees = True
            treenums.append(int(line.split()[-1]))
            if thistree == "":
                continue
            else:
                #temp = open('tmp.txt','w')
                #temp.write(thistree)
                #temp.close()
                trees.append(loadtxt(StringIO(thistree),unpack=True))
                if dopbar:
                    nread = nread + 1
                    pbar.update(nread)
                else:
                    print "Finished tree number {0}".format(treenums[-2])
                thistree = ""
        elif startedtrees:
            thistree = thistree + line
        elif line.startswith('#'):
            comments.append(line.strip('#'))
        else:
            ntrees = int(line)
            nread = 0
            print "Expecting {0} trees in the file".format(ntrees)
            if dopbar:
                widgets = ['Reading trees: ','(',Percentage(),'):',' ',Bar(left='[',right=']',marker='-'),' ',ETA()]
                pbar = ProgressBar(widgets=widgets,maxval=ntrees).start()
            ex = ex + 1
            if ex > 2:
                sys.exit(1)
    if dopbar:
        pbar.finish()
    else:
        print "Reading final tree and saving it..."
    trees.append(loadtxt(StringIO(thistree),unpack=True))
    #print "Removing temporary file tmp.txt"
    #os.remove('tmp.txt')

    assert len(trees) == ntrees

    return [colheads,trees,treenums,comments]


def save_to_hdf5(outname,colheads,trees,treenums,comments):
    '''
    Save tree info to an hdf5 file

    Input:
     outname - hdf5 filename where the trees will be saved
     colheads,trees,treenums,comments - As returned by read_trees()
    '''
    import h5py,os,sys
    print "Saving data to {0}".format(outname)
    f = h5py.File(outname,'a')
    if "RockstarMergerTrees" in f.keys():
        print "File already contains a group named RockstarMergerTrees, so I can't save to it.  Exiting."
        sys.exit(1337)
    t = f.create_group('RockstarMergerTrees')
    for c in comments:
        t.attrs[c] = -1
    for i in range(len(treenums)):
        tre = t.create_group('Tree_{0}'.format(treenums[i]))
        for j in range(len(colheads)):
            if '/' in colheads[j]:
                colheads[j] = colheads[j].replace('/','_')
            tre.create_dataset(colheads[j],data=trees[i][j])
    f.close()


def trees2irate(inname, outname):
    '''
    Read trees in file "inname" and save them to the IRATE file "outname"
    '''
    [colheads,trees,treenums,comments] = read_trees(inname)
    if not outname.endswith('.hdf5'):
        outname = outname + '.hdf5'
    save_to_hdf5(outname,colheads,trees,treenums,comments)

    
## If used as a script ##
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print "usage:  python {0} <inname> <outname>".format(sys.argv[0])
        sys.exit(2)
    inname, outname = sys.argv[1], sys.argv[2] 
    [colheads,trees,treenums,comments] = read_trees(inname)
    if not outname.endswith('.hdf5'):
        outname = outname + '.hdf5'
    save_to_hdf5(outname,colheads,trees,treenums,comments)
