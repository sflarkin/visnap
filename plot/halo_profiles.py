'''Module containing tools for plotting halo radial profiles'''

import os, sys
import h5py
from pylab import *
import visnap
import visnap.plot
from visnap.plot import set_axes
from visnap.functions.mis import find_rpower
#import pdb


def density_profile(halo_object, axes=None, figname='halo_density', fignumber=-1,
                    marker=None, markersize=None, color=None, fontsize=25,
                    showme=1, legendname=None): 
                    
    '''
    Plot the density profile of the given halo

    Input
    
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     marker
     markersize
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
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{log \ Density} \ [\mathrm{M_{\odot}/pc^3}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string, fignumber, fontsize)
    
    # Open halo object
    halo = halo_object
    halo_props = halo.props
    try:
        print 'Plotting density profile of halo %d from %s/%s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot, halo.catalog)
    except AttributeError:
        print 'Plotting density profile of halo %s from %s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot)
    halo.print_properties()
    
    # Generate profiles if not there already
    try:
        halo.profiles
    except AttributeError:
        halo.get_profiles()
    
    rmin = halo.profiles_props['rPower']
    print 'rPower: ', rmin 
    
    # Plot
    h = visnap.h
    r, dens = abs(halo.profiles['r']), halo.profiles['dens']
    simprops = halo.sim_props
    dmName = simprops['dmName']
    if dmName == 'nonCosmo':
        x = r[argwhere(r > rmin)]
        y = dens[argwhere(r > rmin)]/10.0**9
    else:    
        x = r[argwhere(r > rmin)]/h
        y = dens[argwhere(r > rmin)]*h**2/10.0**9

    mpdmS = '%.1e' %  simprops['mpdm']
    mpdmS = mpdmS.split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')
    if color == None: color = 'k'
    if (dmName == 'CDM') or (dmName == 'nonCosmo'):
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            if (dmName == 'nonCosmo'):
                legendname = [str(halo.id)+'\_mpdm'+mpdmS]
            else:    
                legendname = ['halo'+str(int(halo.id))+'\_'
                              +dmName+'\_mpdm'+mpdmS]
    elif dmName == 'WDM':
        if marker == None:  marker = '^'
        if markersize == None: markersize = 10
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_'
                          +simprops['wdmMass']+'kev']
    elif dmName == 'DDM':
        if marker == None:  marker = 's'
        if markersize == None: markersize = 10
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:    
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mp'+mpdmS+'\_'
                          +str(simprops['DDMtau'])+'gyr'+'\_'
                          +str(simprops['DDMvk'])+'km/s' ]
    elif dmName == 'SIDM':
        if marker == None:  marker = '*'
        if markersize == None: markersize = 16
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_s'
                          +str(simprops['sigma_m'])+'\_h'+str(simprops['hsi'])]   
    else: 
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = 'halo'+str(int(halo.id))+'\_' +'\_mpdm'+mpdmS
                            
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


def density_profiles(halo_objects, axes=None, figname='halos_density',fignumber=-1,
                     markers=None, markerssize=None, colors=None, linewidths=None,
                     fontsize=25, showme=1, legendnames=None):
    '''
    Overplot the density profiles of the given halos

    Input:

     halo_objects - A list of halo objects created by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else new lines
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure 

     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()

     markers
     markerssize
     colors      - These will be just passed as given to matplotlib.plot() and
                   they must be lists of equal length or greater than the
                   halo_objects list. If None they will be auto generated 
                                    
     fontsize  - Size of font to use for labels

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendnames - A list of strings with the names to give to the legends. If
                   None they will be auto generated 
     
    Output:

     figure, axes, lines, legends - As returned by matplotlib when creating figures
                                    and lines
     
    '''
    print 'Starting density profiles plot\n'
        
    # Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:    
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{log \ Density} \ [\mathrm{M_{\odot}/pc^3}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)

    # Find number of halos with different zoom id (zoom id given in irate file name)
    simpropss = [halo.sim_props for halo in halo_objects]
    try:
        zoom_ids = [simprops['zoom_id'] for simprops in simpropss]
        NdiffHalos = len(set(zoom_ids))
    except KeyError:
        zoom_ids = arange(len(halo_objects))
        NdiffHalos = zoom_ids.size
       
    # Loop over all halo objects
    lines = []
    legend_names = []
    Nhalos = len(halo_objects)
    for hcount,halo in enumerate(halo_objects):        
        if markers == None: markers = [None]*Nhalos
        if markerssize == None: markerssize = [None]*Nhalos
        if colors is not None:
            color = colors[hcount]
        else:
            color_list = visnap.plot.colors_list
            if NdiffHalos > 1: 
                color = color_list[int(zoom_ids[hcount])] 
            else:
                color = color_list[hcount]

        thisfig, thisax, thisline, thislegend =\
            density_profile(halo, ax, marker=markers[hcount],
                            markersize=markerssize[hcount], 
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


def circular_velocity(halo_object, axes=None, figname='halo_circular_velocity', fignumber=-1,
                    marker=None, markersize=None, color=None, fontsize=25,
                    showme=1, legendname=None): 
                    
    '''
    Plot the circular velocity profile of the given halo

    Input
    
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     marker
     markersize
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
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{Circular Velocity} \ [\mathrm{km/s}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string, fignumber, fontsize)
    
    # Open halo object
    halo = halo_object
    halo_props = halo.props
    try:
        print 'Plotting circular velocity profile of halo %d from %s/%s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot, halo.catalog)
    except AttributeError:
        print 'Plotting circular velocity of halo %s from %s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot)
    halo.print_properties()
    
    # Generate profiles if not there already
    try:
        halo.profiles
    except AttributeError:
        halo.get_profiles()
    
    rmin = halo.profiles_props['rPower']
    print 'rPower: ', rmin 
    
    # Plot
    h = visnap.h
    G = visnap.G
    r, Vc = abs(halo.profiles['r']), sqrt(G*halo.profiles['M_in_r']/abs(halo.profiles['r']))
    simprops = halo.sim_props
    dmName = simprops['dmName']
    if dmName == 'nonCosmo':
        x = r[argwhere(r > rmin)]
        y = Vc[argwhere(r > rmin)]
    else:    
        x = r[argwhere(r > rmin)]/h
        y = Vc[argwhere(r > rmin)]

    mpdmS = '%.1e' %  simprops['mpdm']
    mpdmS = mpdmS.split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')
    if color == None: color = 'k'
    if (dmName == 'CDM') or (dmName == 'nonCosmo'):
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.semilogx(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            if (dmName == 'nonCosmo'):
                legendname = [str(halo.id)+'\_mpdm'+mpdmS]
            else:    
                legendname = ['halo'+str(int(halo.id))+'\_'
                              +dmName+'\_mpdm'+mpdmS]
    elif dmName == 'WDM':
        if marker == None:  marker = '^'
        if markersize == None: markersize = 10
        line = ax.semilogx(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_'
                          +simprops['wdmMass']+'kev']
    elif dmName == 'DDM':
        if marker == None:  marker = 's'
        if markersize == None: markersize = 10
        line = ax.semilogx(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:    
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mp'+mpdmS+'\_'
                          +str(simprops['DDMtau'])+'gyr'+'\_'
                          +str(simprops['DDMvk'])+'km/s' ]
    elif dmName == 'SIDM':
        if marker == None:  marker = '*'
        if markersize == None: markersize = 16
        line = ax.semilogx(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_s'
                          +str(simprops['sigma_m'])+'\_h'+str(simprops['hsi'])]   
    else: 
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.semilogx(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = 'halo'+str(int(halo.id))+'\_' +'\_mpdm'+mpdmS
                            
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


def circular_velocities(halo_objects, axes=None, figname='halos_circular_velocity',fignumber=-1,
                        markers=None, markerssize=None, colors=None, linewidths=None,
                        fontsize=25, showme=1, legendnames=None):
    '''
    Overplot the circular velocity profiles of the given halos

    Input:

     halo_objects - A list of halo objects created by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else new lines
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure 

     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()

     markers
     markerssize
     colors      - These will be just passed as given to matplotlib.plot() and
                   they must be lists of equal length or greater than the
                   halo_objects list. If None they will be auto generated 
                                    
     fontsize  - Size of font to use for labels

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendnames - A list of strings with the names to give to the legends. If
                   None they will be auto generated 
     
    Output:

     figure, axes, lines, legends - As returned by matplotlib when creating figures
                                    and lines
     
    '''
    print 'Starting density profiles plot\n'
        
    # Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:    
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{Circular Velocity} \ [\mathrm{km/s}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)

    # Find number of halos with different zoom id (zoom id given in irate file name)
    simpropss = [halo.sim_props for halo in halo_objects]
    try:
        zoom_ids = [simprops['zoom_id'] for simprops in simpropss]
        NdiffHalos = len(set(zoom_ids))
    except KeyError:
        zoom_ids = arange(len(halo_objects))
        NdiffHalos = zoom_ids.size
       
    # Loop over all halo objects
    lines = []
    legend_names = []
    Nhalos = len(halo_objects)
    for hcount,halo in enumerate(halo_objects):        
        if markers == None: markers = [None]*Nhalos
        if markerssize == None: markerssize = [None]*Nhalos
        if colors is not None:
            color = colors[hcount]
        else:
            color_list = visnap.plot.colors_list
            if NdiffHalos > 1: 
                color = color_list[int(zoom_ids[hcount])] 
            else:
                color = color_list[hcount]

        thisfig, thisax, thisline, thislegend =\
            circular_velocity(halo, ax, marker=markers[hcount],
                              markersize=markerssize[hcount], 
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


def vdisp_profile(halo_object, axes=None, figname='halo_vdisp', fignumber=-1,
                    marker=None, markersize=None, color=None, fontsize=25,
                    showme=1, legendname=None): 
                    
    '''
    Plot the velocity dispersion profile of the given halo

    Input
    
     halo_object - A halo object as returned by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else a new line
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure
               
     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()  
     
     marker
     markersize
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
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{$Velocity \ Dispersion} \ [\mathrm{km/s}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string, fignumber, fontsize)
    
    # Open halo object
    halo = halo_object
    halo_props = halo.props
    try:
        print 'Plotting velocity dispersion profile of halo %d from %s/%s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot, halo.catalog)
    except AttributeError:
        print 'Plotting velocity dispersion profile of halo %s from %s/%s' %\
            (halo.id, halo.irate_file, halo.snapshot)
    halo.print_properties()
    
    # Generate profiles if not there already
    try:
        halo.profiles
    except AttributeError:
        halo.get_profiles()
    
    rmin = halo.profiles_props['rPower']
    print 'rPower: ', rmin 
    
    # Plot
    h = visnap.h
    r, sigv = abs(halo.profiles['r']), halo.profiles['sigv']
    simprops = halo.sim_props
    dmName = simprops['dmName']
    if dmName == 'nonCosmo':
        x = r[argwhere(r > rmin)]
        y = sigv[argwhere(r > rmin)]/10.0**9
    else:    
        x = r[argwhere(r > rmin)]/h
        y = sigv[argwhere(r > rmin)]*h**2/10.0**9

    mpdmS = '%.1e' %  simprops['mpdm']
    mpdmS = mpdmS.split('+')
    mpdmS = mpdmS[0]+mpdmS[1].lstrip('0')
    if color == None: color = 'k'
    if (dmName == 'CDM') or (dmName == 'nonCosmo'):
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            if (dmName == 'nonCosmo'):
                legendname = [str(halo.id)+'\_mpdm'+mpdmS]
            else:    
                legendname = ['halo'+str(int(halo.id))+'\_'
                              +dmName+'\_mpdm'+mpdmS]
    elif dmName == 'WDM':
        if marker == None:  marker = '^'
        if markersize == None: markersize = 10
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_'
                          +simprops['wdmMass']+'kev']
    elif dmName == 'DDM':
        if marker == None:  marker = 's'
        if markersize == None: markersize = 10
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:    
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mp'+mpdmS+'\_'
                          +str(simprops['DDMtau'])+'gyr'+'\_'
                          +str(simprops['DDMvk'])+'km/s' ]
    elif dmName == 'SIDM':
        if marker == None:  marker = '*'
        if markersize == None: markersize = 16
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = ['halo'+str(int(halo.id))+'\_'
                          +dmName+'\_mpdm'+mpdmS+'\_s'
                          +str(simprops['sigma_m'])+'\_h'+str(simprops['hsi'])]   
    else: 
        if marker == None:  marker = '.'
        if markersize == None: markersize = 18
        line = ax.loglog(x, y, color=color, marker=marker, ms=markersize,
                         linewidth=0) 
        if not legendname:
            legendname = 'halo'+str(int(halo.id))+'\_' +'\_mpdm'+mpdmS
                            
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


def vdisp_profiles(halo_objects, axes=None, figname='halos_vdisp',fignumber=-1,
                     markers=None, markerssize=None, colors=None, linewidths=None,
                     fontsize=25, showme=1, legendnames=None):
    '''
    Overplot the velocity dispersion profiles of the given halos

    Input:

     halo_objects - A list of halo objects created by visnap.general.halo.new_halo()

     axes - If None a new figure and axes will be created, else new lines
            will be added to the given axes.

     figname - A string with a name for the created figure file,
               this will be just pased to pylab.savefig(). Set to None if you
               dont want to save the figure 

     fignumber - If the axes are not given a new figure will be creted with
                 this number, if set to -1 it will be given automatically by
                 matplotlib.figure()

     markers
     markerssize
     colors      - These will be just passed as given to matplotlib.plot() and
                   they must be lists of equal length or greater than the
                   halo_objects list. If None they will be auto generated 
                                    
     fontsize  - Size of font to use for labels

     showme - If not 0 the figure will be shown, otherwise it will only be
              saved to a file. The matplotlib figure, axes, and lines
              will be returned for further manipulation 

     legendnames - A list of strings with the names to give to the legends. If
                   None they will be auto generated 
     
    Output:

     figure, axes, lines, legends - As returned by matplotlib when creating figures
                                    and lines
     
    '''
    print 'Starting velocity dispersion profiles plot\n'
        
    # Create plot axes
    if axes:
        fig, ax = axes.figure, axes
    else:    
        xlabel_string = '$\mathrm{r} \ [\mathrm{kpc}]$'
        ylabel_string = '$\mathrm{Velocity \ Dispersion} \ [\mathrm{km/s}]$'
        fig, ax = set_axes.set_standard_axes(xlabel_string, ylabel_string,
                                             fignumber, fontsize)

    # Find number of halos with different zoom id (zoom id given in irate file name)
    simpropss = [halo.sim_props for halo in halo_objects]
    try:
        zoom_ids = [simprops['zoom_id'] for simprops in simpropss]
        NdiffHalos = len(set(zoom_ids))
    except KeyError:
        zoom_ids = arange(len(halo_objects))
        NdiffHalos = zoom_ids.size
       
    # Loop over all halo objects
    lines = []
    legend_names = []
    Nhalos = len(halo_objects)
    for hcount,halo in enumerate(halo_objects):        
        if markers == None: markers = [None]*Nhalos
        if markerssize == None: markerssize = [None]*Nhalos
        if colors is not None:
            color = colors[hcount]
        else:
            color_list = visnap.plot.colors_list
            if NdiffHalos > 1: 
                color = color_list[int(zoom_ids[hcount])] 
            else:
                color = color_list[hcount]

        thisfig, thisax, thisline, thislegend =\
            vdisp_profile(halo, ax, marker=markers[hcount],
                          markersize=markerssize[hcount], 
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
