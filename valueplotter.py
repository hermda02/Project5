import matplotlib.pyplot as plt
import numpy as np
import itertools
import csv

name = 'test_forward.txt'

x=[]
v1 = []

with open(name,'r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in itertools.islice(plots,0,1):
		y = row[0]
		x.append(float(y[2:7]))
		x.append(float(y[8:13]))
		x.append(float(y[14:19]))
		x.append(float(y[20:25]))
		x.append(float(y[26:31]))
		x.append(float(y[32:37]))
		x.append(float(y[38:43]))
		x.append(float(y[44:49]))
		x.append(float(y[50:55]))

with open(name,'r') as csvfile:
	plots = csv.reader(csvfile,delimiter='\t')
	for row in itertools.islice(plots,25,26):
		y = row[0]
		v1.append(float(y[3:14]))
		v1.append(float(y[18:29]))
		v1.append(float(y[33:44]))
		v1.append(float(y[48:59]))
		v1.append(float(y[63:74]))
		v1.append(float(y[78:89]))
		v1.append(float(y[93:104]))
		v1.append(float(y[108:119]))
		v1.append(float(y[123:134]))

name1 = 'test_backward.txt'

v2 = []

with open(name1,'r') as csvfile1:
	plots = csv.reader(csvfile1,delimiter='\t')
	for row in itertools.islice(plots,25,26):
		z = row[0]
		v2.append(float(z[3:14]))
		v2.append(float(z[18:29]))
		v2.append(float(z[33:44]))
		v2.append(float(z[48:59]))
		v2.append(float(z[63:74]))
		v2.append(float(z[78:89]))
		v2.append(float(z[93:104]))
		v2.append(float(z[108:119]))
		v2.append(float(z[123:134]))

name2 = 'test_crank.txt'

v3 = []

with open(name2,'r') as csvfile2:
	plots = csv.reader(csvfile2,delimiter='\t')
	for row in itertools.islice(plots,25,26):
		z = row[0]
		v3.append(float(z[3:14]))
		v3.append(float(z[18:29]))
		v3.append(float(z[33:44]))
		v3.append(float(z[48:59]))
		v3.append(float(z[63:74]))
		v3.append(float(z[78:89]))
		v3.append(float(z[93:104]))
		v3.append(float(z[108:119]))
		v3.append(float(z[123:134]))

plt.title('$\Delta x = 0.1, \Delta t = 0.005$')
plt.xlabel('Rod Length')
plt.ylabel('Temperature')
plt.plot(x,v1,label='FW')
plt.plot(x,v2,label='BW')
plt.plot(x,v3,label='Crank')
plt.legend()
plt.show()
