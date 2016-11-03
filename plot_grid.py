import numpy as np
from plot_custom import *

def cross_slice(data):
     i = 0
     while data['lm_cool'][i]<data['lm_heat_new'][i]:
          i = i+1
     return i    

def add_col1(ax,data,param_file,num,xname,yname1,yname2,yname3,xcoor,ycoor=0.5):
     add_plot(ax[0],param_file,data[xname],np.log10(data[yname1]),xname,yname1,xlabel=False)
     add_plot(ax[1],param_file,data[xname],np.log10(data[yname2]),xname,yname2,xlabel=False)
     add_plot(ax[2],param_file,data[xname],data[yname3],xname,yname3,xlabel=True)
     add_crossline(ax[0],xpoint=data[xname][num],ypoint=data[yname1][num])
     add_crossline(ax[1],xpoint=data[xname][num],ypoint=data[yname2][num])
     add_crossline(ax[2],xpoint=data[xname][num],ypoint=data[yname3][num])
     align_ylabel(ax,xcoor=xcoor,ycoor=ycoor)

def add_col2(ax,data,param_file,num,xname,yname1,yname2,yname3,xcoor,ycoor=0.5):
     add_plot(ax[0],param_file,data[xname],data[yname1],xname,yname1,xlabel=False)
     add_plot(ax[1],param_file,data[xname],data[yname2],xname,yname2,xlabel=False)
     add_plot(ax[2],param_file,data[xname],data[yname3],xname,yname3,xlabel=True)
     add_crossline(ax[0],xpoint=data[xname][num],ypoint=data[yname1][num])
     add_crossline(ax[1],xpoint=data[xname][num],ypoint=data[yname2][num])
     add_crossline(ax[2],xpoint=data[xname][num],ypoint=data[yname3][num])
     align_ylabel(ax,xcoor=xcoor,ycoor=ycoor)

def add_col3(ax,data,param_file,num,xname,yname1,yname2,yname3,xcoor,ycoor=0.5):
     add_plot(ax[1],param_file,data[xname],data[yname2],xname,yname2,xlabel=False)
     add_plot(ax[2],param_file,data[xname],data[yname3],xname,yname3,xlabel=True)
     add_crossline(ax[0],xpoint=data[xname][num],ypoint=data[yname1][num])
     add_crossline(ax[1],xpoint=data[xname][num],ypoint=data[yname2][num])
     add_crossline(ax[2],xpoint=data[xname][num],ypoint=data[yname3][num])
     align_ylabel(ax,xcoor=xcoor,ycoor=ycoor)

def plot_grid(ax,data,param_file):
     num = cross_slice(data)     
     ax1 = [ax[0,1],ax[1,1],ax[2,1]]
     add_col1(ax1,data,param_file,num,'redshift','R_vir','Rcool','m_halo',xcoor=-0.17)
  
     ax2 = [ax[0,0],ax[1,0],ax[2,0]]
     add_col2(ax2,data,param_file,num,'redshift','Hz','lm_bh','logsig1',xcoor=-0.2)

     ax3 = [ax[0,2],ax[1,2],ax[2,2]]
     add_col3(ax3,data,param_file,num,'redshift','lm_cool','lre','m_star',xcoor=-0.15)

     ## plot the mass cooling and heating rate
     add_plot(ax[0,2],param_file,data['redshift'],data['lm_cool'],'redshift','lm_cool',
              xlabel=False,ylim=False,color='blue',label='cooling rate')
     ax[0,2].plot(data['redshift'],data['lm_heat_new'],lw=2,color='red',ls='-',
                  label='heating rate (new)')
     ax[0,2].plot(data['redshift'],data['lm_heat_old'],lw=2,color='red',ls='--',
                  label='heating rate (old)')
     ax[0,2].legend(loc='best',fontsize=12)

     y_min = min(min(data['lm_cool']),min(data['lm_heat_new']),min(data['lm_heat_old']))-0.2
     y_max = max(max(data['lm_cool']),max(data['lm_heat_new']),max(data['lm_heat_old']))+0.2
     ax[0,2].set_ylim(y_min,y_max)
