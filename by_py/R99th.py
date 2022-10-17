# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 15:35:51 2022

计算阈值，说明该事件的极端程度

@author: Lenovo
"""


import os
import pandas as pd
import numpy as np

#%%  读取pre变成dataframe  计算阈值
#! pre 可能有异常值 如何剔除参考 D:\搜狗高速下载\高原极端降水\xiong2019高原极端降水.pdf 

# #读取pre
pre=pd.read_table("F:\\snow_sts_data\\pre1981_2020.txt",sep='\s+',na_values=32700)

# dataframe 中缺失值不影响  np中会影响百分位阈值算法的结果，需去掉
pre.dropna(axis=0, how='any',inplace=True) #删除任何有nan的行 

# 数据选取≥0.1的，因为数据＞0之后最小的就是0.1，所以这里直接筛选＞0的
def fill_existing3(x):
    if x["dailypre"] >0:
        return x["dailypre"]
    else:
        return None
pre.loc[:,'dailypre']= pre.apply(fill_existing3,axis=1)
pre.dropna(axis=0, how='any',inplace=True) #删除任何有nan的行 

# per99=pre['dailypre'].quantile(0.99)    #分位数 不行
#需将数据按照站点、年份分组计算百分位数
pre['year']=pre['time'].astype(str).str[0:4]
pre99_df = pre.groupby(by=['station','year'])['dailypre'].\
    apply(lambda x: np.percentile(x, 99))
pre99_reset=pre99_df.reset_index()
pre99_ave=pre99_reset.groupby(by=['station'])['dailypre'].mean()

#添加位置信息，法二  
path_sta='F:\\snow_sts_data\\after_quality_control\\tp_sta_info_by2014.txt'
sta=pd.read_table(path_sta,sep = ",",usecols=['station','lon','lat'])
sta.set_index('station', inplace=True) # column 改为 index
per99=pd.concat([sta,pre99_ave],axis=1)

per99.dropna(axis=0, how='any',inplace=True) #删除任何有nan的行
a=per99['dailypre'].mean()

# per99_save=per99.reset_index()
# per99_save.to_csv("H:\\d\\data\\precipitation\\per99_1981_2020.txt",index = False,
#                 sep=' ',na_rep=32700)


# %%  读取rashmi日降水  

# rashmi=pd.read_table('H:\\d\\data\\precipitation'+
#                       '\\6h\\sta_r_s_total_py2600_2800.txt',sep="\s+")
# daily=rashmi.loc[:,['tp_sta','prep_oneday']]
# daily.set_index('tp_sta', inplace=True) # column 改为 index
# d=pd.concat([daily,per99],axis=1)
# d.dropna(axis=0, how='any',inplace=True) #删除任何有nan的行 

# def get_label(x):
#     if x['prep_oneday']>x['dailypre']: 
#         return 1
#     else:
#         return 0
# d.loc[:, 'label'] = d.apply(get_label, axis=1) 
# count=d['label'].sum()
# ave=d['dailypre'].mean()


