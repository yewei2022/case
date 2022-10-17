import os
import numpy as np
from netCDF4 import Dataset
dir = 'D:\\ncl_related\\data\\hgt' #路径
ff = os.listdir(dir) #读取目录下的所有文件
print(ff)#输出文件名称
for item in ff:
    nc_file = dir + "\\" + item
    print(nc_file)
# 得到文件的绝对路径，逐个读取文件夹中的nc文件
nc_obj = Dataset(nc_file)#看文件里有什么
print(nc_obj)