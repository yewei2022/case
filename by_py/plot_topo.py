# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 09:50:07 2022

@author: Lenovo
"""

import netCDF4 as nc
import xarray as xr
import numpy as np
fn = xr.open_dataset('H:\\d\\data\\other\\ETOPO2v2c_f4.nc')
z0 = fn['z']
z1=z0.loc[20:50,75:105]
x0 = fn['x']
x1=x0.loc[75:105]
y0 = fn['y']
y1=y0.loc[20:50]
x2 = np.array(x1)
y2 = np.array(y1)
z2 = np.array(z1)
print(x2)
print(y2)
print(z2)
z3=z2.T
print(z3)
z4=np.where(z3 > 0, z3, 0)


#%%
# import matplotlib.pyplot as plt
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter

# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# surf = ax.plot_surface(x2, y2, z4, cmap=cm.coolwarm,
#                         linewidth=0, antialiased=False)
# # Customize the z axis.
# ax.set_zlim(0, 6000)  # z轴的取值范围
# # ax.zaxis.set_major_locator(LinearLocator(10))
# # ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# plt.show()

#%%
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np
from matplotlib import cm

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(x2, y2, z4, rstride=1, cstride=1, cmap=cm.copper, edgecolor='none')
ax.set_title('surface')
# ax.view_init(0, 90) #旋转角度
# 灰色底换成白色
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
# ax.set_zlim(0, 6000)  # z轴的取值范围

pic_dir="H:\\d\\ncl_related\\pictures\\" 
plt.savefig(pic_dir+'topo.ps', dpi=1000, bbox_inches = 'tight')
# 显示图片
plt.show()