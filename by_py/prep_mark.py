#标记高原站点每6h降水量是无0降雨1降雪2
#2006-2900,6h间隔，读取过去6h天气

import os
import numpy as np
#---------------------------读取高原站点-------------------------------
tp_sta =[]
sta_dir='H:\\d\\data\\precipitation\\'
with open(sta_dir+'TP_sta.txt') as sta_file:
    lines=sta_file.readlines()
    for line in lines:
        data=line.split()
        tp_sta.append(data[0]) 
tp_sta.sort()
# print(len(tp_sta)) 97
# exit()

#------------------------根据站点号读取过去天气------------------------------------ 
#过去天气1,2,现在天气分别在第11,12列 ,19             
micaps_dir='H:\\d\\data\\micaps_rain\\need\\'
f_list = os.listdir(micaps_dir)  # 得到文件夹下的所有文件名称
# print(len(f_list)) 12
# list1= list(range(1,len(f_list)+1)) #[1,...,12]
ww=np.zeros((len(tp_sta),12))
i=0
for stad in tp_sta:
    j=0
    for file in f_list:
        outfile=sta_dir+file+".txt"
        filename1=micaps_dir+file
        #print(filename1)
        with open(filename1, 'r', encoding='GB2312') as f1:#有时是'utf-8'
            next(f1)
            next(f1)
            lines = f1.readlines()  # 按行给lines
            for line in lines:
                data=line.split()
                #print(data)
                if data[0] == stad:
                    #print(1)
                    if (int(data[10])==6 or (data[11])==6): 
                        ww[i,j]=1
                    elif (int(data[10])==7 or int(data[11])==7):
                        ww[i,j]=2
                    elif ((int(data[10])==8 and 85<=int(data[18])<=86)
                        or (int(data[11])==8 and 85<=int(data[18])<=86)):
                        ww[i,j]=2
                    elif (int(data[10])==8 or int(data[11])==8):
                        ww[i,j]=1    
                    else:
                        ww[i,j]=0  
        j+=1
    i=i+1
# print(r)                    

# r12=np.zeros(len(sta))
# for i in range(len(sta)):
#     r12[i]=sum(r[i,:])
# print(r12)

outfile=sta_dir+"phase.txt"
with open(outfile, 'w') as fw:
    for i in range(len(tp_sta)):
        fw.write("%d %d %d %d %d %d %d %d %d %d %d %d %d \n"
                 %(int(tp_sta[i]),ww[i,0],ww[i,1],
                   ww[i,2],ww[i,3],ww[i,4],ww[i,5],ww[i,6],
                   ww[i,7],ww[i,8],ww[i,9],ww[i,10],ww[i,11]))