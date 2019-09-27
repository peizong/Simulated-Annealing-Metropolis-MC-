#!/nfs/apps/Compilers/Python/Anaconda/2.7/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

#labels=np.array(["Al","Cr","Fe","Mn","Ti"])
def plot_lro(labels):
  #labels=np.array(["Ti","Fe","Al"]) #np.array(["Fe","Al","Ti"])
  data=np.genfromtxt("lro.dat",skip_header=0)
  for i in range(1,len(data[0])):
   if (i !=2) or (i!=5):
    plt.plot(data[:,0],data[:,i],label=labels[i-1])
  plt.legend()
  plt.xlabel("Temperature/K")
  plt.ylabel("Long range order parameter")
  plt.savefig("lro.png",dpi=300)
  plt.show()
def plot_sro(labels):
  #labels=np.array(["Ti","Fe","Al"]) #np.array(["Fe","Al","Ti"])
  mlabels=[]
  data=np.genfromtxt("sro.dat",skip_header=0)
  for i in range(0,len(labels)):
    a,b=3,2
    idx=str(a)+str(b)+str(i+1)
    plt.subplot(idx)
    for j in range(0,len(labels)):
      #mlabels.append(labels[i]+"-"+labels[j])
      new_label=labels[i]+"-"+labels[j]
      plt.plot(data[:,0],data[:,1+i*len(labels)+j],label=new_label) #mlabels[i*len(labels)+j])
    plt.legend()
    plt.xlabel("Temperature/K")
    plt.ylabel("Short range order parameter")
  plt.savefig("sro.png",dpi=300)
  plt.show()
#labels=np.array(["Ti","Fe","Al"])
labels=np.array(["Al","Co","Cr","Fe","Ni"])
#plot_sro(labels)
plot_lro(labels)
