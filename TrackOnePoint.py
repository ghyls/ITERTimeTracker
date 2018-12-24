#Type here the number corresponding to the FIELD_ID you want to track
#------------------------------------------------------------------
FIELD_ID = 4
smooth_curve = 1

#Set the following variables to "auto" if you want the title to be set 
#automatically. If not, you can know what you are plotting by looking at the 
#table displayed on the terminal.
plot_title = u"Evolution of $B_{z}$"
plot_ylabel = u"$B_z$ ($T$)"
#------------------------------------------------------------------

f = open("trackpoint.txt", 'r')

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# the data inside this file has been generated and has a particulat format. 
# We will first take the data from the file: time and lines:

time = []       #list containing the actual values of time

lines = []      #the lines of the text, splitted by spaces.

for line in f:
    line = list(line.split())
    lines.append(line)

names = [k for k in lines[0][1:]] #first element is "time"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Now we print table with names and values corresponding with the scalars

space = ' '
max_len = max([len(elem) for elem in names])
print("\n\n Available FIELD_ID's:\n") 
print("|   NAME "+(max_len-13)*' '+"   |FIELD_ID |")  
print("|"+max_len*'-'+"+---------|")

for i, elem in enumerate(names):

    sp_units = max_len-len(elem)
    sp = space*sp_units
    print('|'+elem + sp + '     ' + str(i) + '    |')
print("\n\n")
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#we take now the scalars from the file

scalars = [[]] * (len(lines[0])-1)

for line in lines[1:]:
    time.append(float(line[0]))
    for i in range(len(lines[0])-1):
        scalars[i] = scalars[i] + [float(line[i+1])]
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#we plot everything now.

#we will use matplotlib's module "pyplot"
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from numpy import linspace

figure = plt.figure()
ax1 = figure.add_subplot(111)


#we will now interpolate the curve, in order to make the plot much fancier.
resolution = 1000
timenew = linspace(float(time[0]), float(time[-1]), resolution)
ynew = spline(time, scalars[FIELD_ID], timenew)


if not smooth_curve:
    ax1.plot(time, scalars[FIELD_ID], 'b', linewidth=0.5)
else:
    ax1.plot(timenew, ynew, 'b', linewidth=0.5)
    
ax1.scatter(time, scalars[FIELD_ID])

#title and labels
ax1.set_title (names[FIELD_ID] + " vs t", fontsize = 20)
ax1.set_ylabel (names[FIELD_ID], fontsize = 16)

if plot_title != "auto":
    ax1.set_title(plot_title, fontsize=16)
if plot_ylabel != "auto":
    ax1.set_ylabel (plot_ylabel, fontsize = 16)

ax1.set_xlabel (u"time (s)", fontsize = 16)

#number of ticks on the x axis
plt.locator_params(axis='x', nbins=32)

ax1.grid(1,'major', linestyle = '--') 

plt.show()
