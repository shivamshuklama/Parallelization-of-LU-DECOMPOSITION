


import matplotlib
from matplotlib import pyplot as plt

# seq 8.74337   5.47258 3.72565 4.95312  4.61332 5.60211 8.13843

# x axis values
x = [2,4,6,8,16,32]
# corresponding y axis values
y = [5.47258 ,3.72565 ,4.95312  ,4.61332 , 5.60211  , 8.13843]

# plotting the points
plt.plot(x, y)

# naming the x axis
plt.xlabel('NUMBER OF THREADS')
# naming the y axis
plt.ylabel('EXECUTION TIME(SECONDS)')

# giving a title to my graph
plt.title('(taking matrix size 1000)NUMBER OF THREADS vs EXECUTION TIME')


plt.grid(True)

# function to show the plot
#plt.show()


import matplotlib.pyplot as plt1
   
x1 = [2,4,6,8,16,32]
y1 = [0.806, 0.587 , 0.294 , 0.315 , 0.26 , 0.17]

New_Colors = ['green','blue','purple','brown','teal' , 'red']


plot2 = plt.figure(2)
plt1.bar(x1, y1, color=New_Colors)
plt1.title('Pthread (taking matrix size 1000) NUMBER OF THREADS vs PARALLEL EFFICIENCY', fontsize=14)
plt1.xlabel('NUMBER OF THREADS', fontsize=14)
plt1.ylabel('PARALLEL EFFICIENCY', fontsize=14)
plt1.grid(True)

plt.show()
