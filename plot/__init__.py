'''VISNAP modules containing plotting tools'''
 
from pylab import rc, seterr

oldset = seterr(all='ignore') # turn warnings off

#Change some of the default matplotlib configuration options

rc("text.latex",preamble=r'\boldmath')
rc("text",usetex=True)
#rc("mathtext",fontset='custom')
#rc("mathtext",default='bf')
rc("axes", linewidth=3.0)
rc(("xtick.major","ytick.major"),size=15)
rc(("xtick.minor","ytick.minor"),size=8)
rc(("xtick","ytick"),labelsize='x-large') 

#Other definitions useful for the plotting modules

colors_list = ['k','b','g','r','c','m','y','#003300','#990099','#660033','#336666','#663366','#330033','#333333','#666666','#999999' ]
line_styles = ['-','--',':','-.',':',':',':',':',':',':']
marker_styles = ['','','','','o','*','^','+','s']


