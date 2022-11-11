import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F, TFile
from ROOT import gROOT, gBenchmark
import os


#Using a list here instead of a numpy array as I'm pretty sure that lists are faster at appending which we're doing
f = open("data/e2E2e3E3_H250ISR1BS1_LR_addH.lhe_tauDecaysigsig_mixingAngle0.hepevt")

line = f.readline()
evt_num = 0

writer = hep.WriterAscii("data/LR0.hepmc")

while line:
    #print(line)
    numP = int(line)
    px = []
    py = []
    pz = []
    en = []
    m = []
    pid = []
    sta = []
    parents = []


    evt = hep.GenEvent()
     
    for i in range(numP):
        line = f.readline()
        properties = line.split(' ')
        #print(properties)
        fv = hep.FourVector(float(properties[6]),float(properties[7]),float(properties[8]),float(properties[9]))
        px.append(float(properties[6]))
        py.append(float(properties[7]))
        pz.append(float(properties[8]))
        en.append(float(properties[9]))
        m.append(fv.m())
        pid.append(int(properties[1]))
        sta.append(int(properties[0]))
        parents.append([int(properties[2]),int(properties[3])])
    
    #IMPORTANT!! Must be on pyhepmc ~2.7.3 or something. This function is not present in 2.0.0
    evt.from_hepevt(evt_num, px, py, pz, en, m, pid, sta, parents)
    print(evt)
    writer.write_event(evt)    
    line = f.readline()
    evt_num+=1

