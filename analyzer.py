import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F, TFile
from ROOT import gROOT, gBenchmark, Math
import os


#Set to UPDATE to leave previous files there
#Open the .root file in RECREATE mode in order to wipe all histograms previously there
myfile = TFile( 'RootFiles/startup.root', 'RECREATE' )

#
def histogram(directory, name, description):
    hist1 = TH1F( name + "_im", "Invariant Mass of Muon Pairs", 100, 0, 100 )

    directory = os.fsencode(directory)

    reader = None
    cutReader = []

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".hepmc"):
            #needs the .decode, otherwise the directory is in byte format
            reader = hep.ReaderAscii(directory.decode("utf-8") + '/'+ filename)
            cutReader.append(reader)
    print("Found " + str(len(cutReader)) + " files to read.")
    print("Beginning event processing:")

    events = []
    evt = hep.GenEvent()
    endCondition = False
    j = 0
    readerNum = 0
    
    print(dir(cutReader[0]))

    while(not endCondition):
        evt = hep.GenEvent()
        reader.read_event(evt)
        j = j + 1
        if (j % 1000 == 0):
            print(j)
            print(evt)
        if reader.failed():
            readerNum += 1
            #if this is the last reader, end the loop
            if readerNum == len(cutReader):
                endCondition = True
        else:
            events.append(evt)

    print(len(events))
    
    print(events[10].particles[10].parents)
    print(dir(events[10].particles[1]))
    for e in events:
        muon_candidates =[]
        antimuon_candidates = []
        for candidate in e.particles:
            print(candidate.pid)
            if candidate.pid == 13:
                muon_candidates.append(candidate)
            if candidate.pid == -13:
                antimuon_candidates.append(candidate)
        muon = muon_candidates[0]
        antimuon = antimuon_candidates[0]
        for i in muon_candidates:
            if i.momentum.pt() > muon.momentum.pt():
                muon = i
        for i in antimuon_candidates:
            if i.momentum.pt() > antimuon.momentum.pt():
                antimuon = i
        momentum = (muon.momentum + antimuon.momentum)
        hist1.Fill(momentum.m())
        #lorentz = Math.LorentzVector()
        #lorentz.SetPx(momentum.px())
        #lorentz.SetPy(momentum.py())
        #lorentz.SetPz(momentum.pz())
        #lorentz.SetE(momentum.e())

    hist1.Write()
histogram("data","Initial study","H tau tau events")
myfile.Close()
