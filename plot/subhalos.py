'''Module containing tools for the analysis of subhalos'''

import sys
import h5py
from pylab import *
import visnap
from visnap.plot import set_axes
from visnap.general import irate_file_mod
import pdb


def subhalo_function(halo_object, Rmass=-1, vcut_factor=0.05, rcut_factor=1.0,
                     min_npart=100, axes=None, figname='subhalo_function',
                     fignumber=-1, linestyle=None, linewidth=None, color=None,
                     fontsize=25, showme=1, legendname=None): 
    '''
    Plot the subhalo Vmax/Mass function of the given host

    Input
    
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     Rmass - If not -1 the mass function is plotted instead. Rmass can be set as
             'Rvir' or as a numerical value corresponding to the radius in [kpc]
             at which the mass of the subhalos is taken
     
     vcut_factor - Only subhalos with vmax > vcut_factor*host_Vmax will be included

     rcut_factor - Only subhalos within rcut_factor*host_Rvir will be included

     min_npart - Only subhalos with at least min_npart number of particles
                 will be included 

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     linestyle
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
    host = halo_object
    host_props = host.props
    print 'Plotting subhalo function of host halo %d from %s/%s/%s' %\
        (host.id, host.irate_file, host.snapshot, host.catalog)
    host.print_properties()

    # Select subhalos
    subhalos = host.get_subhalos(vcut_factor, rcut_factor, min_npart)
    print 'Got %d subhalos' % len(subhalos)

    # Get some properties
    h = visnap.h
    if Rmass == 'Rvir':
        Mvir_sub = array([sub.props['Mvir'] for sub in subhalos])
        histreturn = ax.hist(Mvir_sub/h,
                             bins=logspace(log10(Mvir_sub.min()/h),
                                           log10(1.1*Mvir_sub.max()/h),num=50),
                             cumulative=-1, histtype='step', log=True, visible=False)
        M_Rmass = Mvir_sub/h
    elif Rmass != -1:
        M_Rmass = -1*ones(len(subhalos))
        r_power = -1*ones(len(subhalos))
        unresolved_subs=0
        below_rpow_subs=0
        for i,sub in enumerate(subhalos): 
            sub.get_profiles()
            rp_sub = sub.profiles['r']
            Mp_sub = sub.profiles['M_in_r'][(abs(rp_sub)/h < Rmass)
                                            & (rp_sub != 0)]
            if Mp_sub.size:
                M_Rmass[i] = Mp_sub.max()
                r_power[i] = sub.profiles_props['rPower'] 
            else: unresolved_subs += 1           
            if r_power[i]/h > Rmass: below_rpow_subs += 1
     
        if (len(subhalos) - unresolved_subs) < 3:
            print "There were only %d subhalos for which the mass profile was "\
                "resolved within %g kpc" % (len(subhalos)-unresolved_subs, Rmass)
            sys.exit()       
        if unresolved_subs > 0:
                print "Warning: There where %d subhalos for which r_power > %g kpc\n" % (bad_subs,Rmass)                
        
        M_Rmass = M_Rmass[M_Rmass != -1]/h
        histreturn = ax.hist(M_Rmass,
                             bins=logspace(log10(M_Rmass.min()),log10(1.1*M_Rmass.max()),num=50), 
                             cumulative=-1, histtype='step', log=True, visible=False)
    else:
        Vmax_sub = array([sub.props['Vmax'] for sub in subhalos])
        histreturn = ax.hist(Vmax_sub,
                             bins=logspace(log10(Vmax_sub.min()),log10(1.1*Vmax_sub.max()),
                                           num=50),
                             cumulative=-1,histtype='step',log=True,visible=False) 
       
    # Plot
    N,bins = histreturn[0], histreturn[1]

    simprops = host.sim_props
    dmName = simprops['dmName']
    mpdmS = '%.1e' %  simprops['mpdm']
    mpdmS = mpdmS.split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')
    if color == None: color = 'k'
    if dmName == 'CDM':
        if linestyle == None: linestyle = 'solid'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:    
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS]
    elif dmName == 'WDM':
        if linestyle == None: linestyle = ':'
        if linewidth == None: linewidth = 5
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_'
                          +simprops['wdmMass']+'kev']
    elif dmName == 'DDM':
        if linestyle == None: linestyle = '-.'
        if linewidth == None: linewidth = 5
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:    
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mp'+mpdmS+'\_'
                          +str(simprops['DDMtau'])+'gyr'+'\_'
                          +str(simprops['DDMvk'])+'km/s' ]
    elif dmName == 'SIDM':
        if linestyle == None: linestyle = '--'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_s'
                          +str(simprops['sigma_m'])+'\_h'+str(simprops['hsi'])]   
    else: 
        if linestyle == None: linestyle = 'solid'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = 'halo'+str(int(host.id))+'\_' +'\_mpdm'+mpdmS
                            
    legend = ax.legend(line,legendname, loc=3,prop=dict(size=0.8*fontsize,),
                       labelspacing=0.1) 
    legend.draw_frame(False)

    for tline in ax.get_xticklines(): tline.set_markeredgewidth(3)
    for tline in ax.get_yticklines(): tline.set_markeredgewidth(3)
    for tline in ax.xaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    for tline in ax.yaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    
    if figname: savefig(figname)
    if showme: show()

    return fig, ax, line[0], legendname[0]


def subhalo_functions(halo_objects, Rmass=-1, vcut_factor=0.05,
                      rcut_factor=1.0,  min_npart=100, axes=None,
                      figname='subhalo_functions', fignumber=-1,
                      linestyles=None, linewidths=None, colors=None, 
                      fontsize=25, showme=1, legendnames=None):
    '''
    Plot the subhalo Vmax/Mass function of the given hosts

    Input
    
     halo_objects - A list of  halo objects created by visnap.general.halo.new_halo()

     Rmass - If not -1 the mass function is plotted instead. Rmass can be set as
             'Rvir' or as a numerical value corresponding to the radius in [kpc]
             at which the mass of the subhalos is taken
     
     vcut_factor - Only subhalos with vmax > vcut_factor*host_Vmax will be included

     rcut_factor - Only subhalos within rcut_factor*host_Rvir will be included

     min_npart - Only subhalos with at least min_npart number of particles
                 will be included 

     axes - If None a new figure and axes will be created, else new lines
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     linestyles
     linewidths
     colors  -   These will be just passed as given to matplotlib.plot() and
                 they must be lists of equal length or greater than the
                 halo_objects list. If None they will be auto generated  

     fontsize  - Size of font to use for labels   

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendnames - A list strings with the name to give to the legends. If
                  None they will auto generated
              
    Output:

     figure, axes, lines, legends - As returned by matplotlib when creating figures
                                    and lines
    ''' 
    print 'Starting subhalo fucntions plot\n'

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

    # Find number of halos with different zoom id (zoom id given in irate file name)
    simpropss = [halo.sim_props for halo in halo_objects]
    zoom_ids = [simprops['zoom_id'] for simprops in simpropss]
    NdiffHalos = len(set(zoom_ids))

    # Loop over all halo objects
    lines = []
    legend_names = []
    Nhalos = len(halo_objects)
    for hcount,halo in enumerate(halo_objects):
        if linestyles == None: linestyles = [None]*Nhalos
        if linewidths == None: linewidths = [None]*Nhalos
        if colors:
            color = colors[hcount]
        else:
            color_list = visnap.plot.colors_list
            if NdiffHalos > 1: 
                color = color_list[int(zoom_ids[hcount])] 
            else:
                color = color_list[hcount]

        thisfig, thisax, thisline, thislegend =\
            subhalo_function(halo, Rmass=Rmass, vcut_factor=vcut_factor,
                             rcut_factor=rcut_factor, min_npart=min_npart,
                             axes=ax, linestyle=linestyles[hcount], 
                             linewidth=linewidths[hcount],
                             color=color, fontsize=fontsize,
                             figname=None, showme=0)       
        
        lines.append(thisline)
        legend_names.append(thislegend)
        print '\n'
        # End of loop over halo objects
                     
    if not legendnames: legendnames = legend_names   
   
    legend = ax.legend(lines, legendnames, loc=3,
                       prop=dict(size=0.8*fontsize,), labelspacing = 0.1)
    legend.draw_frame(False)

    for tline in ax.get_xticklines(): tline.set_markeredgewidth(3)
    for tline in ax.get_yticklines(): tline.set_markeredgewidth(3)
    for tline in ax.xaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    for tline in ax.yaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    
    if figname: savefig(figname)
    if showme: show()

    return fig, ax, lines, legendnames    



def subhalo_2DRfunction(halo_object, projection_axis=(1,0,0),
                        vcut_factor=0.05, rcut_factor=1.0,
                        min_npart=100, x_axis='R',axes=None, figname='subhalo_2DRfunction',
                        fignumber=-1, linestyle=None, linewidth=None, color=None,
                        fontsize=25, showme=1, legendname=None) :
    '''
    Plot the number of subhalos as a function of the 2d distance R from the projection_axis 

    Input:
     
     halo_object - A halo object as returned by visnap.general.halo.new_halo()  
   
     projection_axis - The direction of the projection axis. Can be set as 'major',
                       'intermediate' or 'minor' (principal axes of the 
                       moment of intertia tensor). It can also take a vector in
                       the form (x, y, z).
 
     vcut_factor - Only subhalos with vmax > vcut_factor*host_Vmax will be included

     rcut_factor - Only subhalos within rcut_factor*host_Rvir will be included

     min_npart - Only subhalos with at least min_npart number of particles
                 will be included 

     x_axis - The x_axis can be 'R' or 'Sigma', where Sigma is the average
              enclosed mass surface density             

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     linestyle
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
    
    #Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:
        if x_axis=='Sigma':
            xlabel_string = '$\mathrm{\Sigma}$'
        else:
            xlabel_string = '$\mathrm{R} / \mathrm{R_{Rvir}}$'
        ylabel_string = '$N(\ < \ R)$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)
    
    # Open halo object
    host = halo_object
    host_props = host.props
    print 'Plotting subhalo 2d R function of host halo %d from %s/%s/%s' %\
        (host.id, host.irate_file, host.snapshot, host.catalog)
    host.print_properties()

    # Select subhalos
    subhalos = host.get_subhalos(vcut_factor, rcut_factor, min_npart)
    print 'Got %d subhalos' % len(subhalos)

    # Set normalized vector indicating the projection axes
    if isscalar(projection_axis) and ('Rockstar' in host.catalog):
        print 'At the moment I can only take a vector for the projection_axis '\
            'argument if the halo is from a Rockstar catalog'
        sys.exit()

    if projection_axis == 'major':
        eproj = host.props['Ea']
    elif projection_axis == 'intermediate':
        eproj = host.props['Eb']
    elif projection_axis == 'minor':
        eproj = host.props['Ec']
    else:
        eproj = projection_axis

    eproj = eproj/sqrt(dot(eproj,eproj))  # normalize 
    
    # Find 2d distance from the projection axis of subhalos
    hostCenter = host.props['Center']
    rsub = [sub.props['Center']-hostCenter for sub in subhalos]
    rsub_mag = array([sqrt(dot(r,r)) for r in rsub])
    rsub_dot_eproj = array([dot(r,eproj) for r in rsub])
    Rproj = sqrt(rsub_mag*rsub_mag - rsub_dot_eproj*rsub_dot_eproj)

    if 'Mpc' in host.units['Center']['unitname']:
         Rproj = Rproj*1000

    # Find Sigma if needed
    if x_axis=='Sigma':
        from scipy.interpolate import interp1d
        get2dprofiles = 0
        #get 2-d profiles if not available
        try:
            profiles2d_eproj = host.profiles2d_projection_axis
            if list(profiles2d_eproj) != list(eproj): get2dprofiles = 1
        except AttributeError:
            get2dprofiles = 1
    
        if get2dprofiles:
            host.get_profiles(project_along=eproj, only2d=True)

        Sigma_R = interp1d(host.profiles2d['R'], host.profiles2d['avgSdens'], kind='cubic')
        Sigma_sub = [float(Sigma_R(thisRproj)) for thisRproj in Rproj]

    #plot
    if x_axis=='Sigma':
        histreturn = ax.hist(Sigma_sub,
                             bins=logspace(log10(min(Sigma_sub)),
                                           log10(max(Sigma_sub)),15),
                             cumulative=1, histtype='step', log=True,
                             visible=False)
    else:
        hostRvir = host.props['Rvir']
        histreturn = ax.hist(Rproj/hostRvir,
                             bins=logspace(log10(min(Rproj/hostRvir)), log10(rcut_factor),15),
                             cumulative=1, histtype='step', log=True,
                             visible=False)
                             
    N, bins = histreturn[0], histreturn[1]
    
    simprops = host.sim_props
    dmName = simprops['dmName']
    mpdmS = '%.1e' %  simprops['mpdm']
    mpdmS = mpdmS.split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')
    if color == None: color = 'k'
    if dmName == 'CDM':
        if linestyle == None: linestyle = 'solid'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:    
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS]
    elif dmName == 'WDM':
        if linestyle == None: linestyle = ':'
        if linewidth == None: linewidth = 5
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_'
                          +simprops['wdmMass']+'kev']
    elif dmName == 'DDM':
        if linestyle == None: linestyle = '-.'
        if linewidth == None: linewidth = 5
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:    
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mp'+mpdmS+'\_'
                          +str(simprops['DDMtau'])+'gyr'+'\_'
                          +str(simprops['DDMvk'])+'km/s' ]
    elif dmName == 'SIDM':
        if linestyle == None: linestyle = '--'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = ['halo'+str(int(host.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_s'
                          +str(simprops['sigma_m'])+'\_h'+str(simprops['hsi'])]   
    else: 
        if linestyle == None: linestyle = 'solid'
        if linewidth == None: linewidth = 4
        line = ax.loglog(bins[:-1], N, color=color, lw=linewidth, ls=linestyle)
        if not legendname:
            legendname = 'halo'+str(int(host.id))+'\_' +'\_mpdm'+mpdmS
                            
    legend = ax.legend(line,legendname, loc=3,prop=dict(size=0.8*fontsize,),
                       labelspacing=0.1) 
    legend.draw_frame(False)

    for tline in ax.get_xticklines(): tline.set_markeredgewidth(3)
    for tline in ax.get_yticklines(): tline.set_markeredgewidth(3)
    for tline in ax.xaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    for tline in ax.yaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    
    if figname: savefig(figname)
    if showme: show()

    return fig, ax, line[0], legendname[0]



def subhalo_2DRfunctions(halo_objects, projection_axis=(1,0,0),
                         vcut_factor=0.05, rcut_factor=1.0,
                         min_npart=100, x_axis='R', axes=None, figname='subhalo_2DRfunctions',
                         fignumber=-1, linestyles=None, linewidths=None, colors=None,
                         fontsize=25, showme=1, legendnames=None) :
    '''
    Plot the number of subhalos as a function of the 2d distance R from the projection_axis 

    Input:
     
     halo_object - A halo object as returned by visnap.general.halo.new_halo()  
   
     projection_axis - The direction of the projection axis. Can be set as 'major',
                       'intermediate' or 'minor' (principal axes of the 
                       moment of intertia tensor). It can also take a vector in
                       the form (x, y, z).
 
     vcut_factor - Only subhalos with vmax > vcut_factor*host_Vmax will be included

     rcut_factor - Only subhalos within rcut_factor*host_Rvir will be included

     min_npart - Only subhalos with at least min_npart number of particles
                 will be included 

     x_axis - The x_axis can be 'R' or 'Sigma', where Sigma is the average
             enclosed mass surface density             

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     linestyles
     linewidths
     colors  -  These will be just passed as given to matplotlib.plot() and
                they must be lists of equal length or greater than the
                halo_objects list. If None they will be auto generated  

     fontsize  - Size of font to use for labels   

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The the matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendnames - A list strings with the name to give to the legends. If
                  None they will auto generated
              
    Output:

     figure, axes, line, legend - As returned by matplotlib when creating figures
                                  and lines
                         
    '''

    print 'Starting subhalo 2d R fucntions plot\n'

    #Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:
        if x_axis=='Sigma':
            xlabel_string = '$\mathrm{\Sigma}$'
        else:
            xlabel_string = '$\mathrm{R} / \mathrm{R_{vir}}$'
        ylabel_string = '$N(\ < \ R)$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)
        
    # Find number of halos with different zoom id (zoom id given in irate file name)
    simpropss = [halo.sim_props for halo in halo_objects]
    zoom_ids = [simprops['zoom_id'] for simprops in simpropss]
    NdiffHalos = len(set(zoom_ids))

        # Loop over all halo objects
    lines = []
    legend_names = []
    Nhalos = len(halo_objects)
    for hcount,halo in enumerate(halo_objects):
        if linestyles == None: linestyles = [None]*Nhalos
        if linewidths == None: linewidths = [None]*Nhalos
        if colors:
            color = colors[hcount]
        else:
            color_list = visnap.plot.colors_list
            if NdiffHalos > 1: 
                color = color_list[int(zoom_ids[hcount])] 
            else:
                color = color_list[hcount]

        thisfig, thisax, thisline, thislegend =\
            subhalo_2DRfunction(halo, projection_axis=projection_axis, vcut_factor=vcut_factor,
                                rcut_factor=rcut_factor, min_npart=min_npart,
                                axes=ax, linestyle=linestyles[hcount], 
                                linewidth=linewidths[hcount],
                                color=color, fontsize=fontsize,
                                figname=None, showme=0)       
        
        lines.append(thisline)
        legend_names.append(thislegend)
        print '\n'
        # End of loop over halo objects
                     
    if not legendnames: legendnames = legend_names   
   
    legend = ax.legend(lines, legendnames, loc=3,
                       prop=dict(size=0.8*fontsize,), labelspacing = 0.1)
    legend.draw_frame(False)

    for tline in ax.get_xticklines(): tline.set_markeredgewidth(3)
    for tline in ax.get_yticklines(): tline.set_markeredgewidth(3)
    for tline in ax.xaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    for tline in ax.yaxis.get_minorticklines(): tline.set_markeredgewidth(2)
    
    if figname: savefig(figname)
    if showme: show()

    return fig, ax, lines, legendnames    

    
