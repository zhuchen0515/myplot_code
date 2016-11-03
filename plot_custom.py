from astropy.table import Table
import numpy as np
import  matplotlib.ticker as ticker

## custom axis limit according to the min and max value of the data
def get_lim(data,mticks,tick):
     sec = tick/mticks
     temp = np.amin(data)/sec
     v_min = temp.astype(int)*sec-sec
     temp = np.amax(data)/sec
     v_max = temp.astype(int)*sec+sec
     limit =  [v_min,v_max]
     return limit

def get_med_lim(data,tick):
     median = np.median(data)
     limit = [median-tick,median+tick]
     return limit

## custom axis ticks
def set_axis(ax,xmticks,xtick,ymticks,ytick):
     set_xaxis(ax,xmticks,xtick)
     set_yaxis(ax,ymticks,ytick)

def set_xaxis(ax,xmticks,xtick):
     major_locator = ticker.MultipleLocator(xtick)
     minor_locator = ticker.AutoMinorLocator(xmticks)
     ax.xaxis.set_minor_locator(minor_locator)
     ax.xaxis.set_major_locator(major_locator)

def set_yaxis(ax,ymticks,ytick):
     major_locator = ticker.MultipleLocator(ytick)
     minor_locator = ticker.AutoMinorLocator(ymticks)
     ax.yaxis.set_minor_locator(minor_locator)
     ax.yaxis.set_major_locator(major_locator)

##read the file of axis custom, include axis label, major tick, minor tick
def read_custom(param_file,name):
     list_param = Table.read(param_file,format='ascii.commented_header')
     slice = np.where(list_param['data_name']==name)
     axis_min = list_param['min'][slice]
     axis_max = list_param['max'][slice]
     tick = list_param['tick'][slice]
     mticks = list_param['mticks'][slice]
     axis_label = list_param['label'][slice]
     custom_axis = {'axis_min':axis_min[0],'axis_max':axis_max[0],'tick':tick[0],'mticks':mticks[0],
                    'axis_label':axis_label[0]}
     return custom_axis

def read_custom_cbar(param_file,name):
     list_param = Table.read(param_file,format='ascii.commented_header')
     slice = np.where(list_param['data_name']==name)
     cmin = list_param['cmin'][slice]
     cmax = list_param['cmax'][slice]
     cticks = list_param['cticks'][slice]
     clabel = list_param['label'][slice]
     cmap = list_param['cmap'][slice]
     custom_caxis = {'cmin':cmin[0],'cmax':cmax[0],'cticks':cticks[0],'cmap':cmap[0],
                    'clabel':clabel[0]}
     return custom_caxis
 
##If plot grid plot, align the axis label of each column and row
def align_xlabel(ax,xcoor=0.5,ycoor=-0.17):
     for i in range(np.size(ax)):
          ax[i].get_xaxis().set_label_coords(xcoor,ycoor)

def align_ylabel(ax,xcoor=-0.17,ycoor=0.5):
     for i in range(np.size(ax)):
          ax[i].get_yaxis().set_label_coords(xcoor,ycoor)

## plot cross line
def add_crossline(ax,xpoint,ypoint,lw=1.5,color='black',ls=':',**kwargs):
     ax.axvline(x=xpoint,lw=lw,color=color,ls=ls,**kwargs)
     ax.axhline(y=ypoint,lw=lw,color=color,ls=ls,**kwargs)

#custom the x and y axis
def plot_axis(ax,xdata,ydata,xaxis,yaxis,xlabel=True,ylabel=True,xtick=True,ytick=True,
              xlim=True,ylim=True,xlog=False,ylog=False):
     ##set the ticks of axis
     if xtick == True:
          set_xaxis(ax,xaxis['mticks'],xaxis['tick'])
     if ytick == True:
          set_yaxis(ax,yaxis['mticks'],yaxis['tick'])

     ##set the limit of axis
     if xlim == True:
          limit = get_lim(xdata,xaxis['mticks'],xaxis['tick'])
          ax.set_xlim(limit)
     elif xlim == False:
          ax.set_xlim(xaxis['axis_min'],xaxis['axis_max'])
     elif xlim =='med':
          limit = get_med_lim(xdata,xaxis['tick'])
          ax.set_xlim(limit)

     if ylim == True:
          limit = get_lim(ydata,yaxis['mticks'],yaxis['tick'])
          ax.set_ylim(limit)
     elif ylim == False:
          ax.set_ylim(yaxis['axis_min'],yaxis['axis_max'])
     elif ylim =='med':
          limit = get_med_lim(ydata,yaxis['tick'])
          ax.set_ylim(limit)        

     ##set the axis label 
     if xlabel == True:
          ax.set_xlabel(xaxis['axis_label'])
     elif xlabel == False:
          ax.xaxis.set_ticklabels([])

     if ylabel == True:
          ax.set_ylabel(yaxis['axis_label'])
     elif ylabel == False:
          ax.yaxis.set_ticklabels([])

     if xlog == True:
          ax.set_xscale('log')
          ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
     if ylog == True:
          ax.set_yscale('log')
          ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

def set_axis(ax,param_file,xdata,ydata,xname,yname,xlabel=True,ylabel=True,
               xlim=True,ylim=True,xlog=False,ylog=False,xtick=True,ytick=True):
     xaxis = read_custom(param_file,xname)
     yaxis = read_custom(param_file,yname)
     plot_axis(ax,xdata,ydata,xaxis,yaxis,xlabel=xlabel,ylabel=ylabel,xtick=xtick,ytick=ytick,
              xlim=xlim,ylim=ylim,xlog=xlog,ylog=ylog)

def add_plot(ax,param_file,xdata,ydata,xname,yname,xlabel=True,ylabel=True,
             xlim=True,ylim=True,xlog=False,ylog=False,xtick=True,ytick=True,
             lw=2,color='black',ls='-',**kwargs):
     set_axis(ax,param_file,xdata,ydata,xname=xname,yname=yname,xlabel=xlabel,ylabel=ylabel,
              xtick=xtick,ytick=ytick,xlim=xlim,ylim=ylim,xlog=xlog,ylog=ylog)
     ax.plot(xdata,ydata,lw=lw,color=color,ls=ls,**kwargs)

