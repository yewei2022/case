# -*- coding: utf-8 -*-
"""
Created on Thu Nov 18 10:04:09 2021
重新从网盘下载
@author: Lenovo
"""

"""
=================
Advanced Sounding
=================
Plot a sounding using MetPy with more advanced features.
Beyond just plotting data, this uses calculations from `metpy.calc` to find the lifted
condensation level (LCL) and the profile of a surface-based parcel. The area between the
ambient profile and the parcel profile is colored as well.
"""
# see Metpy and Geocat for examples,cut height
# font default fontname='Helvetica'
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import metpy.calc as mpcalc
from metpy.plots import SkewT
from metpy.units import units


plt.rcParams['font.sans-serif']=['SimHei'] # 用来正常显示中文标签
plt.rcParams['axes.unicode_minus']=False # 用来正常显示负号

###########################################
# Upper air data can be obtained using the siphon package, but for this example we will use  即可以使用siphon包获取探空数据；
# some of MetPy's sample data.

pic_dir='H:/d/py_related/pic/skewt/'
file_dir='H:/d/data/t_logp/08102508-3008/'
filename='2700'
col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed']
df = pd.read_table(file_dir+filename+'.txt', usecols=[0, 1, 2, 3, 4, 5],
                 names=col_names,na_values='9999')
print(df)

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed'), how='all'
               ).reset_index(drop=True)

###########################################
# We will pull the data out of the example dataset into individual variables and
# assign units.

##step2:这一步很重要，将array数据转换为 Unit类型的数据；不转就会出错。
##获取温压湿风信息
p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = df['dewpoint'].values * units.degC
wind_speed = df['speed'].values*2.*units.knots
wind_dir = df['direction'].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)

###########################################
# Create a new figure. The dimensions here give a good aspect ratio.
#创建图形，画探空图，rotation表示探空图中等温线与水平面的夹角。

fig = plt.figure(figsize=(9, 9))
skew = SkewT(fig, rotation = 45)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot.
skew.plot(p, T, 'r',linewidth=3.0)  #环境温度垂直廓形
skew.plot(p, Td, 'g',linewidth=3.0) #环境露点垂直廓线
skew.plot_barbs(p, u, v,length=9,linewidth=1,sizes=dict(spacing=0.15,
                                            height=0.5,width=0.25))  #画出水平风场
skew.ax.set_ylim(750, 150) #设置纵轴(气压)范围
skew.ax.set_xlim(-30, 30)   #设置横轴(温度)范围
# skew.ax.set_aspect(110)#调整长宽比

# Calculate LCL height and plot as black dot. Because `p`'s first value is
# ~1000 mb and its last value is ~250 mb, the `0` index is selected for
# `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
# i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
# should be selected.
# #基于地面气压、温度和露点(以1000hPa充当地面)，计算出抬升凝结高度LCL对应的气压层和温度,并在图上标出
# lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
# skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')

# # Calculate full parcel profile and add to plot as black line
# #画出气块绝热路径
# prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
# skew.plot(p, prof, 'k', linewidth=2)

# # Shade areas of CAPE and CIN
# #填充CAPE和CIN区域
# skew.shade_cin(p, T, prof)
# skew.shade_cape(p, T, prof)

# An example of a slanted line at constant T -- in this case the 0
# isotherm
## Add the relevant special lines 
# Add relevant special lines
# Choose starting temperatures in Kelvin for the dry adiabats
t0 = units.K * np.arange(243.15, 473.15, 10)
skew.plot_dry_adiabats(t0=t0, linestyles='solid', colors='gray', linewidth=1.5)

# Choose temperatures for moist adiabats
t0 = units.K * np.arange(281.15, 306.15, 4)
msa = skew.plot_moist_adiabats(t0=t0,
                               linestyles='solid',
                               colors='lime',
                               linewidths=1.5)

#画出某一条温度 = k的等温线(0)
# skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)

# plot mixing ratios
# Choose mixing ratios
w = np.array([0.001, 0.002, 0.003, 0.005, 0.008, 0.012, 0.020]).reshape(-1, 1)
# Choose the range of pressures that the mixing ratio lines are drawn over
p_levs = units.hPa * np.linspace(1000, 400, 7)
skew.plot_mixing_lines(mixing_ratio=w, pressure=p_levs, colors='lime')

# #画出干绝热线、湿绝热线
# skew.plot_dry_adiabats()
# skew.plot_moist_adiabats()
# skew.plot_mixing_lines()

# Show the plot
plt.grid(True,
         which='major',
         axis='both',
         color='tan',
         linewidth=1.5,
         alpha=0.5)

# font = {'family':'Arial'  #'serif', 
# #         ,'style':'italic'
#         ,'weight':'bold'  # 'normal' 
# #         ,'color':'red'
#         ,'size':20
#         }

# plt.title("Nyingchi"+"  200810"+filename,fontsize=13, pad=15)
plt.xlabel("Temperature (°C)",fontsize=18,fontproperties='Helvetica')#或用font
plt.ylabel("P (hPa)",fontsize=18,fontproperties='Helvetica')
plt.tick_params(labelsize=18)
plt.savefig(pic_dir+filename+'.jpg', dpi=750, bbox_inches = 'tight')
# dpi Adjust sharpness
plt.show()