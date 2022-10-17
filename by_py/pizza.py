# /mnt/d/py_related/test.py
# PYNIO里面可以运行
def make_pizza(size,*toppings):
	"""打印顾客点的所有配料"""
	print("\nMaking a "+str(size)+
		"-inch pizza with the following topings:")
	for topping in toppings:
		print("- "+topping)

