#python /mnt/d/py_related/code/making_pizza.py
# PYNIO里面可以运行
# import pizza
# pizza.make_pizza(16,"wewr","ettg","wtdfk")
import matplotlib.pyplot as plt
input_values=[1,2,3,4,5]
squares=[1,4,9,16,25]
plt.plot(input_values, squares,linewidth=5)
#设置图表标题并给坐标轴加上标签
plt.title("Square Numbers",fontsize=24)
plt.xlabel("Value",fontsize=14)
plt.ylabel("Square of Value",fontsize=14)
#设置刻度标记的大小
plt.tick_params(axis="both",labelsize=14)
plt.savefig('/mnt/d/py_related/pic/h.png') 

#print('\n'.join([''.join([('Welcome!'[(x-y)%8]
#	if((x*0.05)**2+(y*0.1)**2-1)**3-(x*0.05)**2*(y*0.1)**3<=0 else' ')
#	for x in range(-30,30)])for y in range(15,-15,-1)]))