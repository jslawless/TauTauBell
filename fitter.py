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

myfile = uproot.open( 'RootFiles/'+'cp_phase_pi_half'+'.root' )


bellhist = myfile["cp_pi_half_bellInequality;1"].to_numpy()
mockbellhist = myfile["cp_pi_half_mockBellInequality;1"].to_numpy()

xbins = np.delete(np.delete(fix_bins(bellhist[2]),0),-1)
yvalue = fix_bins(bellhist[1])
fig, ax = plt.subplots()

bellfitsx = []
mockbellfitsx = []
bellerror = []
bellfitsy = []
mockbellfitsy = []
mockerror = []

for i in range(0,len(yvalue)):
    xdata = np.delete(np.delete(bellhist[0][i],0),-1)
    popt, pcov = opt.curve_fit(func, xbins, xdata,bounds=([-30,0],[30,50]))
    bellfitsx.append(yvalue[i])
    bellfitsy.append(popt[0])
    bellerror.append(np.sqrt(pcov[0][0]))
    xdata = np.delete(np.delete(mockbellhist[0][i],0),-1)
    popt, pcov = opt.curve_fit(func, xbins, xdata,bounds=([-30,0],[30,50]))
    mockbellfitsx.append(yvalue[i])
    mockbellfitsy.append(popt[0])
    mockerror.append(np.sqrt(pcov[0][0]))

#print(bellfits)
ax.errorbar(bellfitsx,bellfitsy,yerr=bellerror, color = 'b',marker='+',ls='none',capsize=1.3,capthick=1.0)
ax.errorbar(mockbellfitsx,mockbellfitsy,yerr=mockerror, color='r',marker='x',ls='none',capsize=1.3,capthick=1.0)
#ax.plot(xbins, func(xbins, *popt), 'g--', label='fit: a=%5.3f,c=%5.3f' % tuple(popt))

#xbins = np.delete(np.delete(bellhist[2],0),-1)
#hep.histplot(xdata,bins=xbins, ax = ax)
#plt.legend()
ax.text(-9.5, 5.5, "Predicted Bell Effect",
        fontdict={'family': 'arial',
                  'color': 'b',
                  'weight': 'normal',
                  'size': 12,
                  }
        )
ax.text(-9.5, 1.2, "Possible Unknown Cutoff",
        fontdict={'family': 'arial',
                  'color': 'r', 
                  'weight': 'normal',
                  'size': 12,
                  }
        )

ax.axhline(0,color='g',linestyle=':')

ax.set(xlabel=r"Space-Time Separation $\Delta t^2 - \Delta \vec{r}^2$ [mm$^2$]",ylabel="Fitted Amplitude")
plt.show()
