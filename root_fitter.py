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


#def func(x, a, b,c):
#    return a * np.cos(b*x) + c

def func(x, a, c):
    return -a * np.cos(x) + c

myfile = uproot.open( 'RootFiles/'+'cp_phase_0'+'.root' )


bellhist = myfile["cp_0_speedBellInequality;1"].to_numpy()
mockbellhist = myfile["cp_0_speedMockInequality;1"].to_numpy()

xbins = np.delete(np.delete(fix_bins(bellhist[2]),0),-1)
xvalue = fix_bins(bellhist[1])
fig, ax = plt.subplots()

bell_fit_file = open("bell_output.txt")
bell_fit_file.readline()
bell_err_file = open("bell_err.txt")
bell_err_file.readline()
mock_fit_file = open("mock_output.txt")
mock_fit_file.readline()
mock_err_file = open("mock_err.txt")
mock_err_file.readline()

bellfitsx = []
mockfitsx = []
bellerror = []
bellfitsy = []
mockfitsy = []
mockerror = []

print(xvalue)
for i in range(0,len(xvalue)-1):
    bellfitsx.append(xvalue[i])
    mockfitsx.append(xvalue[i])
    bellfitsy.append(float(bell_fit_file.readline()))
    bellerror.append(float(bell_err_file.readline()))
    mockfitsy.append(float(mock_fit_file.readline()))
    mockerror.append(float(mock_err_file.readline()))

ax.errorbar(bellfitsx,bellfitsy,yerr=bellerror, color = 'b',marker='+',ls='none',capsize=1.3,capthick=1.0)
ax.errorbar(mockfitsx,mockfitsy,yerr=mockerror, color='r',marker='x',ls='none',capsize=1.3,capthick=1.0)

print(bellfitsx)

ax.text(0.53, 0.12, "Predicted Bell Effect",
        fontdict={'family': 'arial',
                  'color': 'b',
                  'weight': 'normal',
                  'size': 12,
                  },
        transform=ax.transAxes
        )
ax.text(0.53, 0.6, "v\u03C8  > 2c",
        fontdict={'family': 'arial',
                  'color': 'r', 
                  'weight': 'normal',
                  'size': 12,
                  },
        transform=ax.transAxes
        )

ax.axhline(0,color='g',linestyle=':')

ax.set_xlim([1,5])
ax.set_ylim([-0.05,0.05])
ax.set(xlabel=r"Speed of Mediator $\Delta r / \Delta t$",ylabel="Fitted Amplitude")
plt.show()
