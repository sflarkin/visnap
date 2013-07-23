'''Module containing tools for the analysis of subhalos'''

import h5py
from pylab import *
import visnap.plot
from visnap.plot import set_axes
from visnap.general import translate_filename, find_halo
#import pdb

def zoom_z0_subhalos_incylinder(irate_files, project_along='z', halolist=None,
                                vcut_factor=0.05, rvir_factor=1.0,
                                figname='subhalo_sdensity', fignumber =-1, showme=1,
                                fontsize = 25, min_sub_npart=100, min_npart=100):
    '''
    Overplot the number of subhalos within a cylinder with axis along the
    projection axis defined by project_along

    Input:
     
     irate_files - A list of IRATE files of zoom simulations  
   
     project_along - Set the axis of the cylinder (i.e. line of sight) along
                     this direction, can be x,y, or z

     halolist - If None the halo with the most number particles will be selected as
                the main zoom halo, otherwise a zoomed halo id will be read
                from the catalog name (assuming the hyades convetion) and the main
                zoom halo will be selected by matching the properties of the
                halo with this id in the given halolist (see
                visnap.general.find_halo.find_zoom_halo for more info) 

     vcut_factor - The lowest subhalo vmax to consider in units of the host Vmax

     rvir_factor - Only subhalos within rvir_factor*host_Rvir will be considered

     figname - Name given to the file where the figure is saved (include the
               path on this string to specify the location where to save)

     fignumber - The number given to the created figure, if set to -1
                 it will be given automatically by matplotlib.figure()

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. Note that the matplotlib figure, axes, and lines
              will be returned for further manipulation 

     fontsize - The font size used for labels         

     min_sub_npart - Only subhalos with at leat min_sub_npart particles will be
                     considered

     min_npart - Only halos with min_npart particles will be considered when
                 looking for the main zoomed halo

    Output:

     figure, axes, lines - As returned by matplotlib when creating figures
                           and lines
                         
    '''
    
    print 'Starting zoom subhalo surface density plots\n'

    #get externally defined variables
    h, rhob, rhoCrit = visnap.h, visnap.rho_b, visnap.rho_crit
    colors = visnap.plot.colors_list

    #Create plot axes
    xlabel_string = '$\mathrm{R} / \mathrm{R_{vir}}$'
    ylabel_string = '$N(\ < \ R)$'
    fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string, fignumber, fontsize)

    #Translate file names and set dictionary with sim properties
    nameDicts = [translate_filename.hyades(irate_file) for irate_file in irate_files]
    hids = [nameDict['hid'] for nameDict in nameDicts]
    NdiffHalos = len(set(hids))
    
    #loop over all irate_files
    legendNames = []
    lines = []
    for fcount,irate_file in enumerate(irate_files):
        ahf_flag = 0
        rockstar_flag = 0
        nameDict = nameDicts[fcount] 
        mpdm, hsi, epsilon, sigma, hid = [nameDict['mpdm'], nameDict['hsi'],
                                           nameDict['epsilon'],
                                           nameDict['sigma'], nameDict['hid']]
        data = h5py.File(irate_file,'a')
        values = [value for value in data.values() if "Snapshot" in value.name]
        z0snap = sorted(values)[-1]
        print 'Opening snapshot %s in file %s' % (z0snap.name,irate_file)
        for key in z0snap.keys():
            if 'HaloCatalog' in key:
                catalog_name = key
                break

        if 'AHF' in catalog_name:
            C = z0snap['HaloCatalog_AHF']
            P = C['RadialProfiles']
            ahf_flag = 1
            print "found AHF halo catalog"
        elif 'Rockstar' in catalog_name:
            C = z0snap['HaloCatalog_Rockstar']
            rockstar_flag = 1
            print "found Rockstar halo catalog"
        else:
            print "No AHF or Rockastar halo catalogs were found in this IRATE file"
            sys.exit()
       
        try: 
            hosts = C['Hosts'][...]
        except KeyError:
            try:
                hosts = C['hostHalo'][...]
            except KeyError:
                print "No host/subhalo identifyer found in this IRATE catalog. Now generating"
                routines.subHaloIdentify(C)
                hosts = C['Hosts'][...]

        npart = C['npart'][...]     
        resCut = argwhere(npart > min_npart)[:,0]
        hostCut =  argwhere(hosts == -1)[:,0]
        subCut =  argwhere(hosts != -1)[:,0]
        resCut_sub = argwhere(npart > min_sub_npart)[:,0]
        if ahf_flag:    
            fMhires = C['fMhires'][...]    
            hiresCut = argwhere(fMhires == 1)[:,0]
            finalCut = array([i for i in resCut if (i in hiresCut)])
            finalCut_sub =  array([i for i in resCut_sub if ((i in subCut) & (i in hiresCut))])
            com_offset, mbp_offset, fMhires = [C['com_offset'][...][finalCut],
                                               C['mbp_offset'][...][finalCut],
                                               C['fMhires'][...][finalCut]] 
        else:
            finalCut = array([i for i in resCut])
            finalCut_sub =  array([i for i in resCut_sub if i in subCut])

        Vmax, Rmax, Mvir, Rvir = [C['Vmax'][...][finalCut],
                                  C['Rmax'][...][finalCut],
                                  C['Mvir'][...][finalCut],
                                  C['Rvir'][...][finalCut]] 
        Vmax_sub, Rmax_sub,Mvir_sub, Xc_sub, Yc_sub, Zc_sub  =\
            [C['Vmax'][...][finalCut_sub], C['Rmax'][...][finalCut_sub],
             C['Mvir'][...][finalCut_sub], C['Center'][...][finalCut_sub,0],
             C['Center'][...][finalCut_sub,1], C['Center'][...][finalCut_sub,2]]
                         
        try:
            ID =  C['ID'][...][finalCut]
        except KeyError:
            print 'No halo ID data found in catalog, most probably it is an old AHF catalog. Will set ID = argument'
            ID = arange(npart.size)
            ID = ID[finalCut]
        npart = C['npart'][...][finalCut]
        if halolist:
            hostarg, dist = find_halo.find_zoom_halo(halolist, int(hid), C, finalCut)
        else:
            hostarg = argwhere(npart == npart.max())[0,0]

        halo_center, halo_velocity = C['Center'][...][finalCut][hostarg], C['Velocity'][...][finalCut][hostarg]
        Rvir_host, Vmax_host = Rvir[hostarg], Vmax[hostarg]
        center_units = C['Center'].attrs['unitname']
        print 'Selcting halo %d with %d particles: Mvir = %.2e, Vmax = %g'\
            %(ID[hostarg],npart[hostarg],Mvir[hostarg]/h,Vmax[hostarg]) 
        print 'Center, Velocity: ',halo_center,', ',halo_velocity

        if ahf_flag:
            print "Host halo properties: Vmax = %g km/s, Rmax = %g kpc, Mvir = "\
                "%.2e Msun, Rvir = %g kpc, com_offset = %g, mbp_offset = %g, "\
                "fMhires = %g\n" % (Vmax[hostarg], Rmax[hostarg]/h,
                                    Mvir[hostarg]/h, Rvir[hostarg]/h,
                                    com_offset[hostarg], mbp_offset[hostarg], fMhires[hostarg])    
        else:
            print "Host halo properties: Vmax = %g km/s, Rmax = %g kpc, Mvir = "\
                "%.2e Msun, Rvir = %g kpc \n" % (Vmax[hostarg],Rmax[hostarg]/h,Mvir[hostarg]/h,Rvir[hostarg]/h)
   
        
        R_sub = sqrt((Xc_sub - halo_center[0])*(Xc_sub - halo_center[0]) +
                     (Yc_sub - halo_center[1])*(Yc_sub - halo_center[1]) +
                     (Zc_sub - halo_center[2])*(Zc_sub - halo_center[2])) 

        if project_along == 'z':
            R_project =  sqrt((Xc_sub - halo_center[0])**2 +
                             (Yc_sub - halo_center[1])**2)
        elif project_along == 'y':
            R_project =  sqrt((Xc_sub - halo_center[0])**2 +
                              (Zc_sub - halo_center[2])**2)
        elif project_along == 'x':
            R_project =  sqrt((Yc_sub - halo_center[1])**2 +
                              (Zc_sub - halo_center[2])**2)
        else:
            print "The value of project_along was not recognized, recognized"\
                "options are: 'x', 'y', or 'z'"
                  
        if 'Mpc' in center_units:
            R_sub = R_sub*1000
            R_project = R_project*1000
        arg_sub = argwhere((R_sub < rvir_factor*Rvir_host) & (R_sub > 0) &
                          (Vmax_sub > vcut_factor*Vmax_host))[:,0]
        
        
        histreturn = ax.hist(R_project[arg_sub]/Rvir_host,
                             bins=logspace(log10(min(R_project[arg_sub]/Rvir_host)), log10(rvir_factor),15),
                             cumulative=1, histtype='step', log=True,
                             visible=False)
                             
        N, bins = histreturn[0], histreturn[1]
        icolor = fcount
        if NdiffHalos > 1: icolor = int(hid)

        if nameDict['dmName'] == 'CDM':
            thisline = ax.loglog(bins[:-1], N, color=colors[icolor], lw=4, ls='solid')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']]
        elif nameDict['dmName'] == 'WDM':
            thisline = ax.loglog(bins[:-1], N, color=colors[icolor], lw=5, ls=':')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'
                            +nameDict['wdmMass']+'Kev']
        elif nameDict['dmName'] == 'DDM':
            thisline = ax.loglog(bins[:-1],N,color=colors[icolor],lw=5,ls='-.')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'
                            +str(nameDict['DDMtau'])+'Gyr'+'\_'+str(nameDict['DDMvk'])+'km/s' ]
        else:
            thisline = ax.loglog(bins[:-1],N,color=colors[icolor],lw=4,ls='--')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'+hid
                            +'\_s'+str(sigma)+'\_h'+str(nameDict['hsi'])]
        lines += thisline
        
        data.close()
        #end of loop over irate_files

    for line in ax.get_xticklines(): line.set_markeredgewidth(3)
    for line in ax.get_yticklines(): line.set_markeredgewidth(3)
    for line in ax.xaxis.get_minorticklines(): line.set_markeredgewidth(2)
    for line in ax.yaxis.get_minorticklines(): line.set_markeredgewidth(2)

    savefig(figname)
    if showme:
        legend = ax.legend(lines,legendNames,loc=3,prop=dict(size=0.8*fontsize,),labelspacing = 0.1)
        legend.draw_frame(False)
        show()

    return fig, ax, lines    

def subhalo_function(halo_object, Rmass=-1, axes=None,
                     figname='subhalo_function', fignumber=-1, line=None,
                     linewidth=None, color=None, fontsize=25, showme=1,
                     legendname=None): 
    '''
    Plot the subhalo Vmax/Mass function of the given host

    Input
    
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     Rmass - If not -1 the mass function is plotted instead. Rmass can be set as
             'Rvir' or as a numerical value corresponding to the radius [kpc]
             at which the mass of the subhalo is taken

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     line
     linewidth
     color  -   These will be just passed as given to matplotlib.plot(),
                if None they will be set automatically 

     fontsize  - Size of font to use for labels   

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The the matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendname - As string with the name to give to the legend. If None it
                  will auto generate one
              
    Output:

     figure, axes, line, legend - As returned by matplotlib when creating figures
                                  and lines
    ''' 

    # Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:    
        if Rmass == -1:
            xlabel_string = '$V_\mathrm{max} \ [\mathrm{km/s}]$'
            ylabel_string = '$N(\ > \ V_\mathrm{max})$'
        elif Rmass == 'Rvir':
            xlabel_string = '$M_\mathrm{vir} \ [M_\odot]$'
            ylabel_string = '$N(\ > \ M_\mathrm{vir})$'
        else:
            Rmass = float(Rmass)
            xlabel_string = '$M_\mathrm{' + str(Rmass) + '} \ [M_\odot]$'
            ylabel_string = '$N(\ > \ M_\mathrm{' + str(Rmass) + '})$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)
        
    # Open halo object
    halo = halo_object
    halo_props = halo.props
    print 'Plotting subhalo fucntion of host halo %d from %s/%s/%s' %\
        (halo.id, halo.irate_file, halo.snapshot, halo.catalog)
    halo.print_properties()

    # Plot


    

def zoom_z0_subhalo_functions(irate_files, Rmass=-1, halolist=None, vcut_factor=0.05,
                              rvir_factor=1.0, figname='subhalo_function',
                              fignumber=-1, showme=1, fontsize=25, min_sub_npart=100, min_npart=100):
    '''
    Overplot the subhalo Vmax or Mass functions at z=0 of the zoomed halos
    found in the given IRATE files 
    
    Input: 
    
     irate_files - A list of IRATE files of zoom simulations

     Rmass - If set the mass function is plotted instead. Rmass can be set as
             'Rvir' or as a numerical value corresponding to the radius [kpc]
             at which the mass of the subhalo is taken         
    
     halolist - If None the halo with the most number particles will be selected as
                the main zoom halo, otherwise a zoomed halo id will be read
                from the catalog name (assuming the hyades convetion) and the main
                zoom halo will be selected by matching the properties of the
                halo with this id in the given halolist (see
                visnap.general.find_halo.find_zoom_halo for more info) 

    
     vcut_factor - The lowest subhalo vmax to consider in units of the host Vmax

     rvir_factor - Only subhalos within rvir_factor*host_Rvir will be considered

     figname - Name given to the file where the figure is saved (include the
               path on this string to specify the location where to save)

     fignumber - The number given to the created figure, if set to -1
                 it will be given automatically by matplotlib.figure()

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. Note that the matplotlib figure, axes, and lines
              will be returned for further manipulation             

     fontsize - The font size used for labels

     min_sub_npart - Only subhalos with at leat min_sub_npart particles will be
                     considered

     min_npart - Only halos with min_npart particles will be considered when
                 looking for the main zoomed halo

    Output:

     figure, axes, lines - As returned by matplotlib when creating figures
                           and lines
     
    '''

    print 'Starting subhalo Vmax/Mass function plot routine\n'

    #get externally defined variables
    h, rhob, rhoCrit = visnap.h, visnap.rho_b, visnap.rho_crit
    colors = visnap.plot.colors_list
    
    #Create plot axis
    if Rmass == -1:
        xlabel_string = '$V_\mathrm{max} \ [\mathrm{km/s}]$'
        ylabel_string = '$N(\ > \ V_\mathrm{max})$'
    elif Rmass == 'Rvir':
        xlabel_string = '$M_\mathrm{vir} \ [M_\odot]$'
        ylabel_string = '$N(\ > \ M_\mathrm{vir})$'
    else:
        Rmass = float(Rmass)
        xlabel_string = '$M_\mathrm{' + str(Rmass) + '} \ [M_\odot]$'
        ylabel_string = '$N(\ > \ M_\mathrm{' + str(Rmass) + '})$'
    fig,ax = set_axes.set_standard_axes(xlabel_string, ylabel_string, fignumber, fontsize) # Main figure
    
    if Rmass != -1:
        # Some extra figures that may be useful
        fig100,ax100 = set_axes.set_standard_axes('$V_\mathrm{max} \ [\mathrm{km/s}]$',
                                                  '$R_\mathrm{max} \ [\mathrm{kpc}]$',
                                                  100, fontsize)
        fig101,ax101 = set_axes.set_standard_axes('$V_\mathrm{max} \ [\mathrm{km/s}]$',
                                                  '$R_\mathrm{power} \ [\mathrm{kpc}]$',
                                                  101, fontsize)
        if Rmass == 'Rvir':
            fig102,ax102 = set_axes.set_standard_axes('$V_\mathrm{max} \ [\mathrm{km/s}]$',
                                                      '$M_\mathrm{vir} \ [M_\odot]$',
                                                      102, fontsize)
        else: 
            fig102,ax102 = set_axes.set_standard_axes('$V_\mathrm{max} \ [\mathrm{km/s}]$', 
                                                      '$M_\mathrm{' + str(Rmass) + '} \ [M_\odot]$',
                                                      102, fontsize)
    
    #Translate file names and set dictionary with sim properties
    nameDicts = [translate_filename.hyades(irate_file) for irate_file in irate_files]
    hids = [nameDict['hid'] for nameDict in nameDicts]
    NdiffHalos = len(set(hids))

    #loop over all irate_files
    legendNames = []
    lines = []
    for fcount,irate_file in enumerate(irate_files):
        ahf_flag = 0
        rockstar_flag = 0
        nameDict = nameDicts[fcount] 
        mpdm, hsi, epsilon, sigma, hid = [nameDict['mpdm'], nameDict['hsi'],
                                           nameDict['epsilon'],
                                           nameDict['sigma'], nameDict['hid']]
        data = h5py.File(irate_file,'a')
        values = [value for value in data.values() if "Snapshot" in value.name]
        z0snap = sorted(values)[-1]
        print 'Opening snapshot %s in file %s' % (z0snap.name,irate_file)
        for key in z0snap.keys():
            if 'HaloCatalog' in key:
                catalog_name = key
                break

        if 'AHF' in catalog_name:
            C = z0snap['HaloCatalog_AHF']
            P = C['RadialProfiles']
            ahf_flag = 1
            print "found AHF halo catalog"
        elif 'Rockstar' in catalog_name:
            if (Rmass != -1) & (Rmass != 'Rvir'):
                print 'Only Vmax and Mvir mass functions are supported for Rockstar catalogs\n'
                sys.exit()
            C = z0snap['HaloCatalog_Rockstar']
            rockstar_flag = 1
            print "found Rockstar halo catalog"
        else:
            print "No AHF or Rockastar halo catalogs were found in this IRATE file"
            sys.exit()
       
        try: 
            hosts = C['Hosts'][...]
        except KeyError:
            try:
                hosts = C['hostHalo'][...]
            except KeyError:
                print "No host/subhalo identifyer found in this IRATE catalog. Now generating"
                routines.subHaloIdentify(C)
                hosts = C['Hosts'][...]

        npart = C['npart'][...]     
        resCut = argwhere(npart > min_npart)[:,0]
        hostCut =  argwhere(hosts == -1)[:,0]
        subCut =  argwhere(hosts != -1)[:,0]
        resCut_sub = argwhere(npart > min_sub_npart)[:,0]
        if ahf_flag:    
            fMhires = C['fMhires'][...]    
            hiresCut = argwhere(fMhires == 1)[:,0]
            finalCut = array([i for i in resCut if (i in hiresCut)])
            finalCut_sub =  array([i for i in resCut_sub if ((i in subCut) & (i in hiresCut))])
            com_offset, mbp_offset, fMhires = [C['com_offset'][...][finalCut],
                                               C['mbp_offset'][...][finalCut],
                                               C['fMhires'][...][finalCut]] 
        else:
            finalCut = array([i for i in resCut])
            finalCut_sub =  array([i for i in resCut_sub if i in subCut])

        Vmax, Rmax, Mvir, Rvir = [C['Vmax'][...][finalCut],
                                  C['Rmax'][...][finalCut],
                                  C['Mvir'][...][finalCut],
                                  C['Rvir'][...][finalCut]] 
        Vmax_sub, Rmax_sub,Mvir_sub, Xc_sub, Yc_sub, Zc_sub  =\
            [C['Vmax'][...][finalCut_sub], C['Rmax'][...][finalCut_sub],
             C['Mvir'][...][finalCut_sub], C['Center'][...][finalCut_sub,0],
             C['Center'][...][finalCut_sub,1], C['Center'][...][finalCut_sub,2]]
 
        if (Rmass != -1) and (ahf_flag):
            rp_sub, Mp_sub, Nencp_sub, avgDensp_sub =\
            [P['r'][...][finalCut_sub], P['M_in_r'][...][finalCut_sub],
            P['npart'][...][finalCut_sub], P['ovdens'][...][finalCut_sub]*rhob]

            rmin_sub = array([rp[rp > 0].min() for rp in rp_sub ])
            rp_sub = abs(rp_sub)
                
        try:
            ID =  C['ID'][...][finalCut]
        except KeyError:
            print 'No halo ID data found in catalog, most probably it is an old AHF catalog. Will set ID = argument'
            ID = arange(npart.size)
            ID = ID[finalCut]
        npart = C['npart'][...][finalCut]
        if halolist:
            hostarg, dist = find_halo.find_zoom_halo(halolist, int(hid), C, finalCut)
        else:
            hostarg = argwhere(npart == npart.max())[0,0]

        halo_center, halo_velocity = C['Center'][...][finalCut][hostarg], C['Velocity'][...][finalCut][hostarg]
        Rvir_host, Vmax_host = Rvir[hostarg], Vmax[hostarg]
        center_units = C['Center'].attrs['unitname']
        print 'Selcting halo %d with %d particles: Mvir = %.2e, Vmax = %g'\
            %(ID[hostarg],npart[hostarg],Mvir[hostarg]/h,Vmax[hostarg]) 
        print 'Center, Velocity: ',halo_center,', ',halo_velocity

        if ahf_flag:
            print "Host halo properties: Vmax = %g km/s, Rmax = %g kpc, Mvir = "\
                "%.2e Msun, Rvir = %g kpc, com_offset = %g, mbp_offset = %g, "\
                "fMhires = %g\n" % (Vmax[hostarg], Rmax[hostarg]/h,
                                    Mvir[hostarg]/h, Rvir[hostarg]/h,
                                    com_offset[hostarg], mbp_offset[hostarg], fMhires[hostarg])    
        else:
            print "Host halo properties: Vmax = %g km/s, Rmax = %g kpc, Mvir = "\
                "%.2e Msun, Rvir = %g kpc \n" % (Vmax[hostarg],Rmax[hostarg]/h,Mvir[hostarg]/h,Rvir[hostarg]/h)
   
        R_sub = sqrt((Xc_sub-halo_center[0])*(Xc_sub-halo_center[0]) +
                     (Yc_sub-halo_center[1])*(Yc_sub-halo_center[1]) +
                     (Zc_sub-halo_center[2])*(Zc_sub-halo_center[2])) 

        if 'Mpc' in center_units: R_sub = R_sub*1000
        arg_sub = argwhere((R_sub < rvir_factor*Rvir_host) & (R_sub > 0) &
                           (Vmax_sub > vcut_factor*Vmax_host))[:,0]
        
        if Rmass == 'Rvir':
            histreturn = ax.hist(Mvir_sub[arg_sub]/h,
                                 bins=logspace(log10(Mvir_sub[arg_sub].min()/h),
                                               log10(1.1*Mvir_sub[arg_sub].max()/h),num=50),
                                 cumulative=-1, histtype='step', log=True, visible=False)
            M_Rmass = Mvir_sub[arg_sub]/h
        elif Rmass != -1:
            M_Rmass = -1*ones(len(arg_sub))
            r_power = -1*ones(len(arg_sub))
            bad_subs=0
            for i,thisarg in enumerate(arg_sub):
                rmin = rmin_sub[thisarg]
                if rmin/h > Rmass: bad_subs += 1
                M_Rmass[i] = Mp_sub[thisarg][(abs(rp_sub[thisarg])/h < Rmass) & (rp_sub[thisarg] != 0)].max()
                r_power[i] = rmin
            if bad_subs > 0:
                print "Wargning: There where %d subhalos for which r_power > %g kpc\n" % (bad_subs,Rmass)                
            M_Rmass/=h
            histreturn = ax.hist(M_Rmass,
                                 bins=logspace(log10(M_Rmass.min()),log10(1.1*M_Rmass.max()),num=50), 
                                 cumulative=-1, histtype='step', log=True, visible=False) 
        else:
            histreturn = ax.hist(Vmax_sub[arg_sub],
                                 bins=logspace(log10(Vmax_sub[arg_sub].min()),log10(1.1*Vmax_sub[arg_sub].max()),
                                               num=50),
                                 cumulative=-1,histtype='step',log=True,visible=False)          
        
        N,bins = histreturn[0], histreturn[1]
        icolor = fcount 
        if NdiffHalos > 1: icolor = int(hid)

        if Rmass != -1:   
            ax100.plot(Vmax_sub[arg_sub],Rmax_sub[arg_sub]/h,'.',ms=15,color=colors[icolor])
            if ahf_flag: ax101.plot(Vmax_sub[arg_sub],rmin_sub[arg_sub]/h,'.',ms=15,color=colors[icolor])
            ax102.loglog(Vmax_sub[arg_sub],M_Rmass,'.',ms=15,color=colors[icolor])
       
        if nameDict['dmName'] == 'CDM':
            thisline = ax.loglog(bins[:-1], N, color=colors[icolor], lw=4, ls='solid')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']]
        elif nameDict['dmName'] == 'WDM':
            thisline = ax.loglog(bins[:-1], N, color=colors[icolor], lw=5, ls=':')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'
                            +nameDict['wdmMass']+'Kev']
        elif nameDict['dmName'] == 'DDM':
            thisline = ax.loglog(bins[:-1],N,color=colors[icolor],lw=5,ls='-.')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'
                            +str(nameDict['DDMtau'])+'Gyr'+'\_'+str(nameDict['DDMvk'])+'km/s' ]
        else:
            thisline = ax.loglog(bins[:-1],N,color=colors[icolor],lw=4,ls='--')
            legendNames += [nameDict['dmName']+'\_mp'+nameDict['mpdm']+'\_'+hid
                            +'\_s'+str(sigma)+'\_h'+str(nameDict['hsi'])]
        lines += thisline
        
        data.close()
        # End of loop over irate_files
    
    for line in ax.get_xticklines(): line.set_markeredgewidth(3)
    for line in ax.get_yticklines(): line.set_markeredgewidth(3)
    for line in ax.xaxis.get_minorticklines(): line.set_markeredgewidth(2)
    for line in ax.yaxis.get_minorticklines(): line.set_markeredgewidth(2)
    
    savefig(figname)
    if showme:
        legend = ax.legend(lines,legendNames,loc=3,prop=dict(size=0.8*fontsize,),labelspacing = 0.1)
        legend.draw_frame(False)
        show()

    return fig, ax, lines
