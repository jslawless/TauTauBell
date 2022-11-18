import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F, TH2F,TFile
from ROOT import gROOT, gBenchmark, TLorentzVector
import os


#Set to UPDATE to leave previous files there
#Open the .root file in RECREATE mode in order to wipe all histograms previously there
myfile = TFile( 'RootFiles/2Dhist.root', 'RECREATE' )

#
def histogram(directory, name, description):
    hist1 = TH1F( name + "_im", "Invariant Mass of Muon Pairs", 50, 0, 100 )
    hist2 = TH1F( name + "_tau_im", "Invariant Mass of tau Pairs", 100, 0, 200 )
    hist3 = TH1F( name + "_combined_fv", "Added Magnitude of Tau Four Vectors", 100, -200, 200 )
    hist4 = TH1F( name + "_combined_boosted_fv", "Added Magnitude of Tau Four Vectors Boosted using Muon Frame", 100, -200, 200 )
    hist2D = TH2F( name + "_im_dr", "Invariant Mass of Z vs Tau - AntiTau, no cut", 100, 0, 120, 100, 0, 120 )
    hist2D.GetXaxis().SetTitle("Invariant Mass (GeV)")
    hist2D.GetYaxis().SetTitle("Tau - Antitau (GeV)")

    directory = os.fsencode(directory)

    reader = None
    cutReader = []

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith(".hepmc"):
            #needs the .decode, otherwise the directory is in byte format
            reader = hep.ReaderAscii(directory.decode("utf-8") + '/'+ filename)
            cutReader.append(reader)

    events = []
    evt = hep.GenEvent()
    endCondition = False
    readerNum = 0
    
    while(not endCondition):
        evt = hep.GenEvent()
        reader.read_event(evt)
        if reader.failed():
            readerNum += 1
            #if this is the last reader, end the loop
            if readerNum == len(cutReader):
                endCondition = True
                break
        muon_candidates =[]
        antimuon_candidates = []
        tau_candidates =[]
        antitau_candidates = []
        for candidate in evt.particles:
            if candidate.pid == 13:
                muon_candidates.append(candidate)
            if candidate.pid == -13:
                antimuon_candidates.append(candidate)
            if candidate.pid == 15:
                tau_candidates.append(candidate)
            if candidate.pid == -15:
                antitau_candidates.append(candidate)
        muon = muon_candidates[0]
        antimuon = antimuon_candidates[0]
        tau = tau_candidates[0]
        antitau = antitau_candidates[0]
        for i in muon_candidates:
            if i.momentum.pt() > muon.momentum.pt():
                muon = i
        for i in antimuon_candidates:
            if i.momentum.pt() > antimuon.momentum.pt():
                antimuon = i
        for i in tau_candidates:
            if i.momentum.pt() > tau.momentum.pt():
                tau = i
        for i in antitau_candidates:
            if i.momentum.pt() > antitau.momentum.pt():
                antitau = i
        momentum = (muon.momentum + antimuon.momentum)
        hist1.Fill(momentum.m())
        muon_boost = TLorentzVector()
        muon_boost.SetPx(momentum.px)
        muon_boost.SetPy(momentum.py)
        muon_boost.SetPz(momentum.pz)
        muon_boost.SetE(momentum.e)
        tau_fourvector = TLorentzVector()
        tau_fourvector.SetPx(tau.momentum.px)
        tau_fourvector.SetPy(tau.momentum.py)
        tau_fourvector.SetPz(tau.momentum.pz)
        tau_fourvector.SetE(tau.momentum.e)
        antitau_fourvector = TLorentzVector()
        antitau_fourvector.SetPx(antitau.momentum.px)
        antitau_fourvector.SetPy(antitau.momentum.py)
        antitau_fourvector.SetPz(antitau.momentum.pz)
        antitau_fourvector.SetE(antitau.momentum.e)
        hist2.Fill((tau.momentum+antitau.momentum).m())
        hist3.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())
        tau_fourvector.Boost((91/125)*muon_boost.BoostVector())
        antitau_fourvector.Boost((91/125)*muon_boost.BoostVector())
        hist4.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())
        if((tau_fourvector + antitau_fourvector).Vect().Mag() > 40 and (abs(momentum.m() - 91.19) < 0.5)):
            print()
            print()
            for candidate in evt.particles:
                print(candidate.pid)
                print(candidate.momentum.e)
                print(candidate.status)
        hist2D.Fill(momentum.m(),(tau_fourvector + antitau_fourvector).Vect().Mag())
    hist1.Write()
    hist2.Write()
    hist3.Write()
    hist4.Write()
    hist2D.Write()
histogram("data","Initial study","H tau tau events")
myfile.Close()
