import numpy as np
import sys
sys.path.append('/data1/zhuchen/SIP/package')
from data_candels.join import join_table

field = ['gds','uds','gdn','cos','egs']
data = join_table(field)
slice = np.where((data['hmag_obs']<=25.5) &\
                 (data['mass']>=9.0) &\
                 (data['sig1']>0) &\
                 (data['galfit_flag']==0) &\
                 (data['redshift']>=0.5) & (data['redshift']<2.5) &\
                 (data['PhotFlag']==0) & \
                 (data['CLASS_STAR']<0.9) &\
                 (data['sf_flag']==1) & (data['delssfr']>=-0.4))

data_obs = data[slice]

redshift = data_obs['redshift']
lmass = data_obs['mass']
sig1 = data_obs['sig1']
lre = data_obs['sma']
