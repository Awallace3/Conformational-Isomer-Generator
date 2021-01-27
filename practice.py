import numpy as np
import random
from mpl_toolkits import mplot3d
from PIL import Image
import matplotlib.pyplot as plt
print('What is the backbone file name?' + ' please add .dat at the end')
backbone = input()
with open('carts.txt','r') as f:
    line = []
    amountoflines = 0
    for i in f.readlines():
        amountoflines += 1
        line.append(i)
        a_array = []
        for i in range(amountoflines):
            a_array.append(line[i].split())
            numberofatoms = a_array[0]
    #    print(update[1:])
    #    print(a_array[1:])
        a_array= a_array[1:]

      #  a_array = np.array(update)
       # print(a_array)
     #   print(a_array[0][1])
        el = []
        x = []
        y = []
        z = []
        for i in a_array:
            el.append(float(i[0]))
            x.append(float(i[1]))
            y.append(float(i[2]))
            z.append(float(i[3]))
        nx = np.array(x)
        ny = np.array(y)
        nz = np.array(z)
        d = np.column_stack((x,y,z,el))
#print(d)




#a.dat, d.dat, benzene,dat
#meth, et, d
#a = ["Spears", "Adele", "NDubz", "Nicole", "Cristina"]
#b = [1, 2, 3, 4, 5]
#c = ["h","k","i"]
electronacceptor = meth
backbone = d
electrondonor = et
electronacceptorfloat = meth.astype(float)
backbonefloat = d.astype(float)
electrondonor = et.astype(float)

fig1= plt.figure(figsize=(12,10))
ax = fig1.add_subplot(1, 2, 2, projection='3d')
for num,i in enumerate(backbonefloat[:,3]):
    #print(i)
    if i == 1.0:
        ax.plot(backbonefloat[num,0],backbonefloat[num,1],backbonefloat[num,2],'ob',linestyle='None')
    elif i == 6.0:
        ax.plot(backbonefloat[num,0],backbonefloat[num,1],backbonefloat[num,2],'ok',linestyle='None')  
#print(plt.show())
#fig2= plt.figure(figsize=(12,10))
ay = fig1.add_subplot(1, 2, 1, projection='3d')
for num,i in enumerate(electronacceptorfloat[:,3]):
   # print(i)
    if i == 1.0:
        ay.plot(electronacceptorfloat[num,0],electronacceptorfloat[num,1],electronacceptorfloat[num,2],'ob',linestyle='None')
    elif i == 6.0:
        ay.plot(electronacceptorfloat[num,0],electronacceptorfloat[num,1],electronacceptorfloat[num,2],'ok',linestyle='None')  



#print(plt.show())

solar = [electronacceptor] + [backbone] + [electrondonor]
#list(solar)
#print(solar[2])
random.shuffle(solar)

#solar = np.concatenate(
#combined = list(zip(a, b, c))

#random.shuffle(combined)

#a[:], b[:], c[:]  = zip(*combined)
#print(solar)










