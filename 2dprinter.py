import matplotlib.pyplot as plt
import numpy as np
import itertools
import csv

name = '2dim_exp_10x10.txt'

data=[]

with open(name,'r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		a = str(row[0])
		b = a.split()
		data.append(b)

print(data)

