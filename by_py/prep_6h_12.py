
"""
读取所有出现过的站点 6h降水量 保存为rain12_py.txt 但是好像也还是97个站点啊 高原站点数
重新从网盘下载

"""

import os
import numpy as np
import pandas as pd

#%% 循环读取数据

# #挑出3日内所有所有文件中出现过的站号，站点，经纬度，6h降水量
# #计算累计降水量,保存为rain12_py.txt
# sta = []

# path='D:/data/precipitation/6h/need_6h/'
# f_list = os.listdir(path)  # 得到文件夹下的所有文件名称
# #print(len(f_list))

# #1 获取所有出现过的站号
# for file in f_list:
#     filename1=path+file
#     #print(filename1)
#     with open(filename1, 'r', encoding='UTF-8') as f1:
#         next(f1)
#         lines = f1.readlines()  # 按行给lines
#         for line in lines:
#             data=line.split()
#             #print(data)
#             if data[4] not in sta:
#                 sta.append(data[4])
                
# sta.sort()
# print(sta)
# print(len(sta))

# #2 根据站号读取相应的变量            
# lat=np.zeros(len(sta))
# lon=np.zeros(len(sta))
# r=np.zeros((len(sta),12))
# i=0
# for stad in sta:
#     j=0
#     for file in f_list:
#         filename1=path+file
#         #print(filename1)
#         with open(filename1, 'r', encoding='UTF-8') as f1:
#             next(f1)
#             lines = f1.readlines()  # 按行给lines
#             for line in lines:
#                 data=line.split()
#                 #print(data)
#                 if data[4] == stad:
#                     #print(1)
#                     lat[i]=float(data[5])
#                     lon[i]=float(data[6])
#                     if float(data[8])!=-1.0:
#                         r[i,j]=float(data[8])
#         j+=1
#     i=i+1
# print(r)                    

# r12=np.zeros(len(sta))
# for i in range(len(sta)):
#     r12[i]=sum(r[i,:])
# print(r12)

# #3 写入文件
# path1='D:/data/precipitation/6h/rain12_py.txt'
# with open(path1, 'w') as fw:
#     for i in range(len(sta)):
#         fw.write("%s %7.2f %7.2f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f \n"
#                  %(sta[i],lat[i],lon[i],r[i,0],r[i,1],
#                    r[i,2],r[i,3],r[i,4],r[i,5],r[i,6],
#                    r[i,7],r[i,8],r[i,9],r[i,10],r[i,11],r12[i]))


#%% 再次读取，重新处理

# info=pd.read_table("E:\\d\\data\\precipitation\\6h\\rain12_py.txt",sep="\s+",
#                    header=None)
# info.columns=['sta','lat','lon','2606','2612','2618','2700','2706','2712',
#               '2718','2800','2806','2812','2818','2900','2600_2900']
# keys=['2606','2612','2618','2700','2706','2712','2718','2800']
# info.loc[:,'2600_2800']=info[keys].sum(1)
# info.to_csv(r"E:\\d\\data\\precipitation\\6h\\rain12_py.txt",index = False,
#                 sep=' ',na_rep=0)
