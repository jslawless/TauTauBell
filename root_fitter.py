import os
import uproot
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import scipy.optimize as opt
import matplotlib.patheffects as path_effects

def fix_bins(bins):
    arr = []
    for i in range(0, len(bins)-1):
        arr.append(float(round((bins[i+1]-bins[i])/2+bins[i],4)))
    return np.array(arr)

def func(x, a, c):
    return -a * np.cos(x) + c

class Plot:
    def __init__(self,na,cl,ma,x,y,desc):
        self.name = na
        self.fitsx = []
        self.error = []
        self.fitsy = []
        self.color = cl
        self.marker = ma
        self.fit_file = open("data/"+na+"_output.txt")
        self.fit_file.readline()
        self.err_file = open("data/"+na+"_err.txt")
        self.err_file.readline()
        self.xpos = x
        self.ypos = y
        self.description = desc

myfile = uproot.open( 'RootFiles/'+'cp_phase_0'+'.root' )

plots = {
    "mock2": Plot("mock2",'g','o',0.53, 0.7, "v\u03C8  > 3c"),
    "bell": Plot("bell",'b','+',0.53, 0.12, "Predicted Bell Effect"),
    "mock": Plot("mock",'r','x',0.53, 0.6, "v\u03C8  > 2c")
        }

smeared_plots = {
    "smmock2": Plot("smmock2",'g','o',0.53, 0.7, "v\u03C8  > 3c, then smeared"),
    "smbell": Plot("smbell",'b','+',0.53, 0.12, "Bell Effect Smeared"),
    "smmock": Plot("smmock",'r','x',0.53, 0.6, "v\u03C8  > 2c, then smeared")
    }

bellhist = myfile["cp_0_speedBellInequality;1"].to_numpy()

xbins = np.delete(np.delete(fix_bins(bellhist[2]),0),-1)
xvalue = fix_bins(bellhist[1])
fig, ax = plt.subplots(2)

for i in range(0,len(xvalue)-1):
    for pl in plots:
        plots[pl].fitsx.append(xvalue[i])
        plots[pl].fitsy.append(float(plots[pl].fit_file.readline()))
        plots[pl].error.append(float(plots[pl].err_file.readline()))
    for pl in smeared_plots:
        smeared_plots[pl].fitsx.append(xvalue[i])
        smeared_plots[pl].fitsy.append(float(smeared_plots[pl].fit_file.readline()))
        smeared_plots[pl].error.append(float(smeared_plots[pl].err_file.readline()))


for pl in plots:
    ax[0].errorbar(plots[pl].fitsx,plots[pl].fitsy,yerr=plots[pl].error, color = plots[pl].color,marker=plots[pl].marker,ls='none',capsize=1.3,capthick=1.0)
    ax[0].text(plots[pl].xpos, plots[pl].ypos, plots[pl].description,
            fontdict={'family': 'arial',
                  'color': plots[pl].color,
                  'weight': 'normal',
                  'size': 12,
                  },
            transform=ax[0].transAxes
            )
for pl in smeared_plots:
    ax[1].errorbar(smeared_plots[pl].fitsx,smeared_plots[pl].fitsy,yerr=smeared_plots[pl].error, color = smeared_plots[pl].color,marker=smeared_plots[pl].marker,ls='none',capsize=1.3,capthick=1.0)
    ax[1].text(smeared_plots[pl].xpos, smeared_plots[pl].ypos, smeared_plots[pl].description,
            fontdict={'family': 'arial',
                  'color': smeared_plots[pl].color,
                  'weight': 'normal',
                  'size': 12,
                  },
            transform=ax[1].transAxes
            )

ax[0].axhline(0,color='g',linestyle=':')

#ax.set_xlim([1,5])
ax[0].set_ylim([-0.05,0.05])
ax[0].set(xlabel=r"Speed of Mediator $\Delta r / \Delta t$",ylabel="Fitted Amplitude")
ax[1].axhline(0,color='g',linestyle=':')

#ax.set_xlim([1,5])
ax[1].set_ylim([-0.05,0.05])
ax[1].set(xlabel=r"Speed of Mediator $\Delta r / \Delta t$",ylabel="Fitted Amplitude")
plt.tight_layout()
plt.show()
