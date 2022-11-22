import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F, TFile
from ROOT import gROOT, gBenchmark
import os


f = open("mixing_angle_data/e2E2e3E3_H250ISR1BS1_LR_addH.lhe_tauDecaysigsig_mixingAngle0.5.hepevt")

line = f.readline()
evt_num = 0

writer = hep.WriterAscii("mixing_angle_data/LR0.5.hepmc")

while line:
#while evt_num < 4: 
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
    children = []
    vx = []
    vy = []
    vz = []
    vt = []

    evt = hep.GenEvent()
     
    for i in range(numP):
        line = f.readline()
        properties = line.split(' ')
        #print(properties)
        px.append(float(properties[6]))
        py.append(float(properties[7]))
        pz.append(float(properties[8]))
        en.append(float(properties[9]))
        m.append(float(properties[10]))
        pid.append(int(properties[1]))
        sta.append(int(properties[0]))
        parents.append([int(properties[2]),int(properties[3])])
        children.append([int(properties[4]),int(properties[5])])
        vx.append(float(properties[11]))
        vy.append(float(properties[12]))
        vz.append(float(properties[13]))
        vt.append(float(properties[14]))

    #IMPORTANT!! Must be on pyhepmc ~2.7.3 or something. This function is not present in 2.0.0
    evt.from_hepevt(evt_num, px, py, pz, en, m, pid, sta, parents, children, vx, vy, vz, vt)
    print(evt)
    writer.write_event(evt)    
    line = f.readline()
    evt_num+=1
print(evt_num)
