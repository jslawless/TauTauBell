import numpy as np

def func(x, a,b):
    return (b) / (a * np.cos(x) +b)

def find_bin(val,bins):
    #these bins represent bin edges. So whichever
    #bin is the first one that its less than, that's the bin that it fits inuuu
    for i in range(0,len(bins)):
        if val < float(bins[i]):
            return i-1
    return (len(bins)-2)

def reweight(angle,speed, amp, bins,offset):
    bin_num = find_bin(speed,bins)
    value = func(angle,float(amp[bin_num]),float(offset[bin_num]))
    return value 
