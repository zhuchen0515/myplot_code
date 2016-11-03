import numpy as np
import sys
from astropy.table import Table

cosmos = sys.argv[1]

def get_name(file_name,num):
     data = Table.read(file_name,format='ascii.commented_header')
     return data['name'][num]

def get_data(dir_name,cosmos,mass):
     if dir_name == 'lin':
          dir_path = '/data1/zhuchen/SIP/project/black_hole/data_lin/out_data_'+cosmos+'/'
     elif dir_name == 'moster':
          dir_path = '/data1/zhuchen/SIP/project/black_hole/data_moster/out_data_'+cosmos+'/'
     elif dir_name == 'dong':
          dir_path = '/data1/zhuchen/SIP/project/black_hole/data_dong/out_data_'+cosmos+'/'

     data_param = '/data1/zhuchen/SIP/project/black_hole/param/out_data.param'
     param = Table.read(data_param, format='ascii.commented_header')
     slice = np.where((param['cosmos']==cosmos) & (param['data']==dir_name) &\
                      (param['mass']==mass))
     n = param[slice]['a']
     data_name = get_name(dir_path+'ini_list.cat',n)
     data_file = dir_path+'out_data/'+data_name[0]+'.cat'
     data = Table.read(data_file, format='ascii.commented_header')
     return data

data_lin_low = get_data('lin',cosmos,10.0)
data_lin_mid = get_data('lin',cosmos,10.5)
data_lin_high = get_data('lin',cosmos,11.0)

data_moster_low = get_data('moster',cosmos,10.0)
data_moster_mid = get_data('moster',cosmos,10.5)
data_moster_high = get_data('moster',cosmos,11.0)

data_dong_low = get_data('dong',cosmos,10.0)
data_dong_mid = get_data('dong',cosmos,10.5)
data_dong_high = get_data('dong',cosmos,11.0)

def com_data(name):
     data_low = [data_lin_low[name],data_moster_low[name],data_dong_low[name]]
     data_mid = [data_lin_mid[name],data_moster_mid[name],data_dong_mid[name]]     
     data_high = [data_lin_high[name],data_moster_high[name],data_dong_high[name]]
     return data_low, data_mid, data_high

