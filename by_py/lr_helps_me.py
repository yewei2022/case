# from liran, read data
import os

filename="E:/data/wcm/505stid.txt"

station = []
with open(filename, 'r', encoding='UTF-8') as f:
    next(f)
    lines = f.readlines()  # 按行给lines
    for line in lines:
        mm = line.split()
        station.append(int(mm[0]))#append代表在数组station和后面连起来


path="E:/data/wcm/PRE/"
f_list = os.listdir(path)  # 得到文件夹下的所有文件名称
#print(len(f_list))
for sta in station:
    stid = []
    year = []
    month = []
    day = []
    rain = []
    for file in f_list:
        filename1=path+file
        #print(filename1)
        with open(filename1, 'r', encoding='UTF-8') as f1:
            lines = f1.readlines()  # 按行给lines
            for line in lines:
                mm = line.split()
                if len(mm) >0 :
                    if int(mm[0]) == sta:
                        stid.append(int(mm[0]))
                        year.append(int(mm[4]))
                        month.append(int(mm[5]))
                        day.append(int(mm[6]))
                        rain.append(int(mm[9]))

    path1="E:/data/wcm/result/"
    outfile= path1+('%5d.dat'%(sta))
    with open(outfile, 'w') as fw:
        for i in range(len(stid)):
            fw.write("%6d %5d %3d %3d %7d\n"%(stid[i],year[i],month[i],day[i],rain[i]))

