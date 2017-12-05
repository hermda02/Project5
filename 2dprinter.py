import matplotlib.pyplot as plt
import numpy as np
import itertools
import csv

name = 'two_dim_exp.txt'

data=[]

with open(name,'r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in plots:
		y = row
		new = y.split()
		print(new)