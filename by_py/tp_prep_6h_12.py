#读取高原站号，经纬度，3日内所有6h降水量
#计算累计降水量，保存为tp_rain12_py.txt
#当时是在Ubuntu 的 PYNIO里运行的 当然这里也可以直接运行
# python /mnt/d/py_related/code/tp_prep_6h_12.py

import os
import numpy as np

#1 读取高原站点
tp_sta =[]
sta_dir='/mnt/d/data/precipitation/6h/TP_sta.txt'
with open(sta_dir) as sta_file:
    lines=sta_file.readlines()
    for line in lines:
        data=line.split()
        tp_sta.append(data[0]) 
tp_sta.sort()

#2 根据站号读取相应的变量  
path='/mnt/d/data/precipitation/6h/need_6h/'
f_list = os.listdir(path)  # 得到文件夹下的所有文件名称
#print(len(f_list))          
lat=np.zeros(len(tp_sta))
lon=np.zeros(len(tp_sta))
r=np.zeros((len(tp_sta),12))
i=0
for stad in tp_sta:
    j=0
    for file in f_list:
        filename1=path+file
        #print(filename1)
        with open(filename1, 'r', encoding='UTF-8') as f1:
            next(f1)
            lines = f1.readlines()  # 按行给lines
            for line in lines:
                data=line.split()
                #print(data)
                if data[4] == stad:
                    #print(1)
                    lat[i]=float(data[5])
                    lon[i]=float(data[6])
                    if float(data[8])!=-1.0:
                        r[i,j]=float(data[8])
        j+=1
    i=i+1

r12=np.zeros(len(tp_sta))
for i in range(len(tp_sta)):
    r12[i]=sum(r[i,:])

#写入文件
path1='/mnt/d/data/precipitation/6h/tp_rain12_py.txt'
with open(path1, 'w') as fw:
    for i in range(len(tp_sta)):
        fw.write("%s %7.2f %7.2f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f \n"
                 %(tp_sta[i],lat[i],lon[i],r[i,0],r[i,1],
                   r[i,2],r[i,3],r[i,4],r[i,5],r[i,6],
                   r[i,7],r[i,8],r[i,9],r[i,10],r[i,11],r12[i]))
