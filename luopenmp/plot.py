



from matplotlib import pyplot as plt

#seq 11.0287 sec 6.2489 4.00389 3.49932 5.47976 5.62361   6.38336 


# x axis values
x = [2,4,6,8,16,32]
# corresponding y axis values s/p
y = [6.2489 ,4.00389 ,3.49932 ,5.47976, 5.62361  , 6.38336]


plot1 = plt.figure(1)

# plotting the points
plt.plot(x, y)

# naming the x axis
plt.xlabel('NUMBER OF THREADS')
# naming the y axis
plt.ylabel('Time Taken For Execution(seconds)')

# giving a title to my graph
plt.title('(taking matrix size 1000) NUMBER OF THREADS vs Time Taken For Execution ')
plt.grid(True)

# function to show the plot
#plt.show()



import matplotlib.pyplot as plt1
   
x1 = [2,4,6,8,16,32]
y1 = [0.8887 , 0.6887 , 0.5252 , 0.335 , 0.327 , 0.288]

New_Colors = ['green','blue','purple','brown','teal' , 'red']


plot2 = plt.figure(2)
plt1.bar(x1, y1, color=New_Colors)
plt1.title(' Openmp (taking matrix size 1000) NUMBER OF THREADS vs PARALLEL EFFICIENCY', fontsize=14)
plt1.xlabel('NUMBER OF THREADS', fontsize=14)
plt1.ylabel('PARALLEL EFFICIENCY', fontsize=14)
plt1.grid(True)

plt.show()
