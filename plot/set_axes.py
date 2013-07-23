'''Module containing tools for setting up plot axes'''

from pylab import figure

def set_standard_axes(xlabel_string, ylabel_string, fig_number=-1,
                      font_size=25, fig_size=(10,8)):
    '''
    Create a figure and set the axes with the corresponding labels

    Input:
     xlabel_string 
     ylabel_string    
     fig_number - The number given to the created figure, if set to -1
                  it will be given automatically by pylab.figure()
     font_size -  The font size used for labels
     fig_size  -  A tuple with the desired figure dimensions

    Output:
     fig - The figure
     ax - The axes
    '''
    if fig_number == -1:
        fig = figure(figsize = fig_size)
    else: 
        fig = figure(fig_number, figsize = fig_size)
    ax = fig.add_subplot(111)
    fig.subplots_adjust(bottom = 0.13,top=0.93)
    ax.set_xlabel(xlabel_string, size=font_size)
    ax.set_ylabel(ylabel_string, size=font_size)
    
    return fig, ax

