import matplotlib.pyplot as plt
import numpy as np
#import os
#os.system('main.exe')
def Func_Plot(n, axes, m):
    str1 = ["P", "rho", "u", "S", "T"]
    str2 = ["result_p.txt", "result_rho.txt", "result_u.txt", "result_S.txt", "result_T.txt"]
    for j in range(m):
        f = open(str2[j], 'r')
        try:
            data_plot = []
            lines = f.readlines()
            #print(len(lines))
            for line in lines:
                data_plot.append([float(x) for x in line.split(' ') if x != "\n"])
        finally:
            f.close()
    # ----------------------------------------------
        #print(len(data_plot[0]), len(data_plot[1]))
        x = data_plot[0]
        t = []
        data_plot = data_plot[1:len(data_plot)]
        for data in data_plot:
            t.append(data[0])
            data.pop(0)
       # print(t)
       # print(len(t))
    
        h = 70
       # print(len(data_plot[h*2]), len(x))
        #print(data_plot)
        for i in range(n):
        #    print(i*h)
        #    print(len(x), len(data_plot[i*h]))
            axes[j, i].plot(x, data_plot[i*h])
            #axes[i].set_title(, fontsize = 10, fontweight = 'heavy', color = 'red')
            axes[j, i].text(-1.2, 0.7, "t="+str(round(t[i*h],2))+" sec",fontsize = 8, fontweight = 'heavy', color = 'orange')
            axes[j, i].set_xlim(x[0]-0.2, x[len(x)-1]+0.2)
            if j==4:
                axes[j, i].set_ylim(0.59, 1.41)
            else:
                axes[j, i].set_ylim(-0.01, 1.5)
                
            axes[j, i].plot([-1, -1],[2, 0], "black", linewidth = 3)
            axes[j, i].plot([1, 1],[2, 0], "black", linewidth = 3)
            axes[j, i].grid()
            if(i == 0):
                axes[j, i].set_ylabel(str1[j],fontsize=15, fontweight='bold')
        y = [0, 0.2, 0.4, 0.6, 0.8, 1]#[x/2 for x in range(0,3)]
        plt.setp(axes[0:3], yticks = y)
        plt.setp(axes[4], yticks = [0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8])
    

#plt.plot([0, 1],[0, 1./0.621059], linewidth = 1)
#plt.plot([0, -1],[0, 1./1.18322], linewidth = 1)
#plt.plot([0, 1],[0, 1./0.68313], linewidth = 1)
#plt.plot([-1, 1],[0.5, 0.5], linewidth = 1)
#plt.legend(["V = 0.621059",'a4 = 1.18322','a1 = 0.68313'])
#plt.show()
n_time = 2
fig, axes = plt.subplots(6,n_time)
Func_Plot(n_time, axes, 5)
f = open("x_t_diagram.txt", 'r')
try:
    data_plot = []
    lines = f.readlines()
    #print(len(lines))
    for line in lines:
        data_plot.append([float(x) for x in line.split(' ') if x != "\n"])
finally:
    f.close()    
#print(data_plot)
speeds = np.array(data_plot[0])
x_plot = np.array(data_plot[1])
#print(speeds[0], x_plot)
t_plot1 =  x_plot/speeds[0]
t_plot2 =  x_plot/speeds[1]
axes[5,0].plot(x_plot, t_plot1, x_plot, t_plot2)
axes[5,0].legend(["Sound speed:"+str(speeds[0]), "Shock wave speed:"+str(speeds[1])])
plt.show()

