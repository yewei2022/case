"""
计算过程累计、日累计的 降雪量，降雨量，降水量
保存为sta_r_s_total_py2600_2800.txt
"""

import os
import pandas as pd
import numpy as np

#%% 法一 循环读取 #根据高原站点的降水相态标记 计算2600-2900、2612-2712的累计降雪量，降雨量，降水量
# # python /mnt/d/py_related/code/rain_snow.py  linux
# # python d:/py_related/code/rain_snow.py   windows
# # 可直接打开Ubuntu运行，里面的python是搭载在anaconda上的，比较方便

# #---------------------------读取高原站点号和降水相态-------------------------------
# sta_dir='/mnt/d/data/precipitation/6h/'
# tp_sta=[]
# with open(sta_dir+'phase.txt') as sta_file:
#     lines=sta_file.readlines()
#     phase=np.zeros((len(lines),12))
#     i=0
#     for line in lines:
#         data=line.split()
#         tp_sta.append(data[0])
#         for j in range(12):
#             phase[i,j]=int(data[j+1])
#         i+=1
# # print(type(tp_sta[0]))
# # exit()

# #------------------------根据站号读取相应的变量------------------------------------               
# rain_file='/mnt/d/data/precipitation/6h/rain12_py.txt'
# r=np.zeros((len(tp_sta),12))
# prep_total=np.zeros(len(tp_sta))
# lat=np.zeros(len(tp_sta))
# lon=np.zeros(len(tp_sta))
# i=0
# for stad in tp_sta:
#     #如果下面这行不放在循环里面而在外面，那么文件只能读取一次
#     with open(rain_file, 'r', encoding='UTF-8') as fr: 
#         lines = fr.readlines()  # 按行给lines
#         for line in lines:
#             data=line.split()
#             if data[0] == stad:
#                 lat[i]=float(data[1])
#                 lon[i]=float(data[2])
#                 prep_total[i]=float(data[15])
#                 for j in range(12):
#                     r[i,j]=float(data[j+3])
#         i=i+1
# # print(r[36,:]) 
# # exit()
# r_total=np.zeros(len(tp_sta))
# s_total=np.zeros(len(tp_sta))
# r_oneday=np.zeros(len(tp_sta))
# s_oneday=np.zeros(len(tp_sta))
# prep_oneday=np.zeros(len(tp_sta))
# for i in range(len(tp_sta)):
#     s12=np.zeros(12)
#     s12=np.where(phase[i,:]==2,r[i,:],0)
#     s_total[i]=sum(s12[:])
#     s_oneday[i]=sum(s12[2:6])
#     prep_oneday[i]=sum(r[i,2:6]) #注意第3-第6，4个数的索引为2:6
#     r_total[i]=prep_total[i]-s_total[i]
#     r_oneday[i]=prep_oneday[i]-s_oneday[i]

# outfile= '/mnt/d/data/precipitation/6h/sta_r_s_total_py.txt'
# with open(outfile, 'w') as fw:
#     for i in range(len(tp_sta)):
#         fw.write("%s %7.2f %7.2f %5.1f %5.1f %5.1f %5.1f %5.1f %5.1f \n"
#                  %(tp_sta[i],lat[i],lon[i],r_total[i],s_total[i],prep_total[i],
#                     r_oneday[i],s_oneday[i],prep_oneday[i]))
        
        
#%% 法二 pandas读取 #根据高原站点的降水相态标记 计算2600-2800、2612-2712的累计降雪量，降雨量，降水量

# phase1=pd.read_table("H:\\d\\data\\precipitation\\6h\\phase.txt",sep="\s+",
#                     header=None)
# phase1.columns=['sta','2606','2612','2618','2700','2706','2712',
#               '2718','2800','2806','2812','2818','2900']
# tp_sta=phase1.sta.tolist()

# prep=pd.read_table("H:\\d\\data\\precipitation\\6h\\rain12_py.txt",sep="\s+")
# prep1=prep[prep.sta.isin(phase1.sta)].reset_index(drop=True)

# sta_info=prep1[['sta','lat','lon']]
# phase1.drop(columns = ['sta','2806','2812','2818','2900'],inplace = True)
# prep1.drop(columns = ['sta','lat','lon','2600_2900','2600_2800','2806','2812',
#                       '2818','2900'],inplace = True)

# cols=prep1.shape[1]
# r=np.array(prep1)
# phase=np.array(phase1)


# r_total=np.zeros(len(tp_sta))
# s_total=np.zeros(len(tp_sta))
# prep_total=np.zeros(len(tp_sta))
# r_oneday=np.zeros(len(tp_sta))
# s_oneday=np.zeros(len(tp_sta))
# prep_oneday=np.zeros(len(tp_sta))

# for i in range(len(tp_sta)):
#     s12=np.zeros(cols) 
#     prep_total[i]=sum(r[i,:])
#     s12=np.where(phase[i,:]==2,r[i,:],0)
#     s_total[i]=sum(s12[:])
#     s_oneday[i]=sum(s12[2:6])
#     prep_oneday[i]=sum(r[i,2:6]) #索引[2:6]区间表示索引为2-5 即第3-6列这4列数据
#     r_total[i]=prep_total[i]-s_total[i]
#     r_oneday[i]=prep_oneday[i]-s_oneday[i]
# names=['tp_sta','r_total','s_total','prep_total','r_oneday','s_oneday','prep_oneday']
# df=pd.DataFrame(list(zip(tp_sta,r_total,s_total,prep_total,r_oneday,s_oneday,
#                           prep_oneday)),columns=names)

# # 添加位置信息 法一  循环添加
# lat_lon=[]
# for i in tp_sta:
#     ii=int(i)
#     lat_lon.append(sta_info.loc[sta_info["sta"] == ii,['lat','lon']])
# lat_lon_df = pd.concat(lat_lon,ignore_index=True)
# df1=pd.concat([lat_lon_df,df],axis=1)

# # df1.to_csv("H:\\d\\data\\precipitation\\6h\\sta_r_s_total_py2600_2800.txt",
# #             index = False,sep=' ',columns=['tp_sta','lat','lon','r_total','s_total',
# #                                           'prep_total','r_oneday',
# #                                           's_oneday','prep_oneday'])


#%% 测试 降雪站点数有多少

# phase1=pd.read_table("D:\\case\\data\\precipitation\\6h\\sta_r_s_total_py2600_2800.txt",
#                       sep="\s+")
# a1=phase1[phase1['s_total']>=0.1]
# a2=phase1[phase1['s_oneday']>=0.1]
# print('累计降雪站次: {}, 单日降雪站次: {}'.format(len(a1),len(a2)))



#%% pandas 计算2606-2718

# phase1=pd.read_table("H:\\d\\data\\precipitation\\6h\\phase.txt",sep="\s+",
#                     header=None)
# phase1.columns=['sta','2606','2612','2618','2700','2706','2712',
#               '2718','2800','2806','2812','2818','2900']
# tp_sta=phase1.sta.tolist()

# prep=pd.read_table("H:\\d\\data\\precipitation\\6h\\rain12_py.txt",sep="\s+")
# prep1=prep[prep.sta.isin(phase1.sta)].reset_index(drop=True)

# #开始处理
# sta_info=prep1[['sta','lat','lon']]
# phase1.drop(columns = ['sta','2606','2800','2806','2812','2818','2900'],inplace = True)
# prep1.drop(columns = ['sta','lat','lon','2600_2900','2600_2800','2606',
#                       '2800','2806','2812','2818','2900'],inplace = True)

# cols=prep1.shape[1]
# r=np.array(prep1)
# phase=np.array(phase1)
# # print(r)

# r_total=np.zeros(len(tp_sta))
# s_total=np.zeros(len(tp_sta))
# prep_total=np.zeros(len(tp_sta))

# r_oneday=np.zeros(len(tp_sta))
# s_oneday=np.zeros(len(tp_sta))
# prep_oneday=np.zeros(len(tp_sta))

# for i in range(len(tp_sta)):
#     s12=np.zeros(cols)
#     prep_total[i]=sum(r[i,:])
#     s12=np.where(phase[i,:]==2,r[i,:],0)
#     s_total[i]=sum(s12[:])
#     s_oneday[i]=sum(s12[1:5]) #索引1-5区间表示索引为1-4 即第2-5列这4列数据
#     prep_oneday[i]=sum(r[i,1:5]) 
#     r_total[i]=prep_total[i]-s_total[i]
#     r_oneday[i]=prep_oneday[i]-s_oneday[i]
# names=['tp_sta','r_total','s_total','prep_total','r_oneday','s_oneday','prep_oneday']
# df=pd.DataFrame(list(zip(tp_sta,r_total,s_total,prep_total,r_oneday,s_oneday,
#                           prep_oneday)),columns=names)

# # 添加位置信息 法一  循环添加
# lat_lon=[]
# for i in tp_sta:
#     ii=int(i)
#     lat_lon.append(sta_info.loc[sta_info["sta"] == ii,['lat','lon']])
# lat_lon_df = pd.concat(lat_lon,ignore_index=True)
# df1=pd.concat([lat_lon_df,df],axis=1)
# df1.to_csv("H:\\d\\data\\precipitation\\6h\\sta_r_s_total_py2606_2718.txt",
#             index = False,sep=' ',columns=['tp_sta','lat','lon','r_total','s_total',
#                                           'prep_total','r_oneday',
#                                           's_oneday','prep_oneday'])


#%% 将55690 56227 56434 三站2600-2800的6h降水量储存

# prep=pd.read_table("H:\\d\\data\\precipitation\\6h\\rain12_py.txt",sep="\s+")
# prep.drop(columns = ['lat','lon','2600_2900','2600_2800',
#                       '2806','2812','2818','2900'],inplace = True)

# Cona=prep[prep['sta']==55690]
# Bome=prep[prep['sta']==56227]
# Zayu=prep[prep['sta']==56434]
# df=pd.concat([Cona,Bome,Zayu])
# df_T = pd.DataFrame(df.values.T,columns=df.index,index=df.columns) #转置
# df_T.drop('sta',axis=0,inplace=True)
# df_T.to_csv("H:\\d\\data\\precipitation\\6h\\prep_6h_3sta.txt",
#             index = False,sep=' ',header=None)



