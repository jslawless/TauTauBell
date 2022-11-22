import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F,TH2F, TFile
from ROOT import gROOT, gBenchmark, TLorentzVector
import os


#Set to UPDATE to leave previous files there
#Open the .root file in RECREATE mode in order to wipe all histograms previously there
myfile = TFile( 'RootFiles/tau_cross_piplus_mixing_angle.root', 'RECREATE' )

#
def histogram(directory, name, description):
    hist1 = TH1F( name + "_im", "Invariant Mass of Muon Pairs", 50, 0, 100 )
    hist2 = TH1F( name + "_tau_im", "Invariant Mass of tau Pairs", 100, 0, 200 )
    hist3 = TH1F( name + "_combined_fv", "Added Magnitude of Tau Four Vectors", 100, -200, 200 )
    hist4 = TH1F( name + "_combined_boosted_fv", "Added Magnitude of Tau Four Vectors Boosted using Muon Frame", 100, -200, 200 )
    hist5 = TH1F( name + "_angleBetweenUnitVecs", "Angle between Unit vectors of plane, tau cross pi plus, with a non SM mixing angle", 50, -5, 5) 
    hist6 = TH1F( name + "_angleBetweenTaus", "Angle between Tau Threevectors", 50, -5,5)
    hist2D = TH2F( name + "_im_dr", "Invariant Mass of Z vs Tau - AntiTau, no cut", 100,81, 101, 100, 0, 12 )
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
        piplus_candidates = []
        piminus_candidates = []
        for candidate in evt.particles:
            if candidate.pid == 13:
                muon_candidates.append(candidate)
            if candidate.pid == -13:
                antimuon_candidates.append(candidate)
            if candidate.pid == 15:
                tau_candidates.append(candidate)
            if candidate.pid == -15:
                antitau_candidates.append(candidate)
            if candidate.pid == 211:
                piplus_candidates.append(candidate)
            if candidate.pid == -211:
                piminus_candidates.append(candidate)
        muon = muon_candidates[0]
        antimuon = antimuon_candidates[0]
        tau = tau_candidates[0]
        antitau = antitau_candidates[0]
        piplus = piplus_candidates[0]
        piminus = piminus_candidates[0]
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
        for i in piplus_candidates:
            if i.momentum.pt() > piplus.momentum.pt():
                piplus = i
        for i in piminus_candidates:
            if i.momentum.pt() > piminus.momentum.pt():
                piminus = i
        momentum = (muon.momentum + antimuon.momentum)
        hist1.Fill(momentum.m())
        muon_boost = TLorentzVector()
        muon_boost.SetPtEtaPhiM(momentum.pt(),momentum.eta(),momentum.phi(),125)
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
        piplus_fourvector = TLorentzVector()
        piplus_fourvector.SetPx(piplus.momentum.px)
        piplus_fourvector.SetPy(piplus.momentum.py)
        piplus_fourvector.SetPz(piplus.momentum.pz)
        piplus_fourvector.SetE(piplus.momentum.e)
        piminus_fourvector = TLorentzVector()
        piminus_fourvector.SetPx(piminus.momentum.px)
        piminus_fourvector.SetPy(piminus.momentum.py)
        piminus_fourvector.SetPz(piminus.momentum.pz)
        piminus_fourvector.SetE(piminus.momentum.e)
        hist2.Fill((tau.momentum+antitau.momentum).m())
        hist3.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())
        tau_fourvector.Boost(muon_boost.BoostVector())
        antitau_fourvector.Boost(muon_boost.BoostVector())
        piplus_fourvector.Boost(muon_boost.BoostVector())
        piminus_fourvector.Boost(muon_boost.BoostVector())
        hist4.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())
        #u1 = tau_fourvector.Vect().Cross(piminus_fourvector.Vect()).Unit()
        #u2 = antitau_fourvector.Vect().Cross(piplus_fourvector.Vect()).Unit()
        u1 = tau_fourvector.Vect().Cross(piplus_fourvector.Vect()).Unit()
        u2 = antitau_fourvector.Vect().Cross(piminus_fourvector.Vect()).Unit()
        angleforSign = u1.Cross(u2).Unit()
        sign = -1
        if (angleforSign.Angle(tau_fourvector.Vect()) < (math.pi/2) ):
            sign = 1
        hist5.Fill(sign*u1.Angle(u2))
        hist6.Fill(tau_fourvector.Vect().Angle(antitau_fourvector.Vect()))
        hist2D.Fill(momentum.m(),(tau_fourvector + antitau_fourvector).Vect().Mag())
    myfile.Write()
histogram("mixing_angle_data","Initial study","H tau tau events")
myfile.Close()
