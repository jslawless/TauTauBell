import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TFile, TArrow, TPaveLabel
from ROOT import gROOT, gBenchmark, TLorentzVector
from ROOT import kMagenta, kTRUE, kFALSE
import os


#Set to UPDATE to leave previous files there
#Open the .root file in RECREATE mode in order to wipe all histograms previously there
myfile = TFile( 'RootFiles/event_display.root', 'RECREATE' )


def final_state_momentum_display():

    c1 = TCanvas("c1","Final State Momentum Display", 0, 0, 200, 200)
    reader = hep.ReaderAscii('./data/single_event.hepmc')

    events = []
    evt = hep.GenEvent()
    endCondition = False
    reader.read_event(evt)
    a1 = TArrow()
    label = TPaveLabel()
    muon = None
    antiMuon = None
    piplus = None
    piminus = None
    for particle in evt.particles: 
        if particle.status == 1 and particle.pid != 22:
            if particle.pid == 13:
                muon = particle
            if particle.pid == -13:
                antimuon = particle
            if particle.pid == 211:
                piplus = particle
            if particle.pid == -211:
                piminus = particle 
            #added 0.5 to the particle's end position to account for the center being at 0.5
            #not sure best way to do that
            a1.DrawArrow(0.5,0.5,(particle.momentum.px/200)+0.5,(particle.momentum.py/200)+0.5,0.007,"|>") 
            label.DrawPaveLabel((particle.momentum.px/200)+0.51,(particle.momentum.py/200)+0.51,(particle.momentum.px/200)+0.53,(particle.momentum.py/200)+0.53,str(particle.pid))

    c1.Write()

    momentum = (muon.momentum + antimuon.momentum)
    muon_boost = TLorentzVector()
    muon_boost.SetPtEtaPhiM(momentum.pt(),momentum.eta(),momentum.phi(),125)
    muon_fourvector = TLorentzVector()
    muon_fourvector.SetPx(muon.momentum.px)
    muon_fourvector.SetPy(muon.momentum.py)
    muon_fourvector.SetPz(muon.momentum.pz)
    muon_fourvector.SetE(muon.momentum.e)
    antimuon_fourvector = TLorentzVector()
    antimuon_fourvector.SetPx(antimuon.momentum.px)
    antimuon_fourvector.SetPy(antimuon.momentum.py)
    antimuon_fourvector.SetPz(antimuon.momentum.pz)
    antimuon_fourvector.SetE(antimuon.momentum.e)
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

    muon_fourvector.Boost(muon_boost.BoostVector())
    antimuon_fourvector.Boost(muon_boost.BoostVector())
    piplus_fourvector.Boost(muon_boost.BoostVector())
    piminus_fourvector.Boost(muon_boost.BoostVector())
    c2 = TCanvas("c2","Boosted Final State",0,0,200,200)

    a1.DrawArrow(0.5,0.5,(muon_fourvector.Px()/200)+0.5,(muon_fourvector.Py()/200)+0.5,0.007,"|>")
    label.DrawPaveLabel((muon_fourvector.Px()/200)+0.51,(muon_fourvector.Py()/200)+0.51,(muon_fourvector.Px()/200)+0.54,(muon_fourvector.Py()/200)+0.54,"13")

    a1.DrawArrow(0.5,0.5,(antimuon_fourvector.Px()/200)+0.5,(antimuon_fourvector.Py()/200)+0.5,0.007,"|>")
    label.DrawPaveLabel((antimuon_fourvector.Px()/200)+0.51,(antimuon_fourvector.Py()/200)+0.51,(antimuon_fourvector.Px()/200)+0.54,(antimuon_fourvector.Py()/200)+0.54,"-13")
    a1.DrawArrow(0.5,0.5,(piplus_fourvector.Px()/200)+0.5,(piplus_fourvector.Py()/200)+0.5,0.007,"|>")
    label.DrawPaveLabel((piplus_fourvector.Px()/200)+0.51,(piplus_fourvector.Py()/200)+0.51,(piplus_fourvector.Px()/200)+0.54,(piplus_fourvector.Py()/200)+0.54,"211")
    a1.DrawArrow(0.5,0.5,(piminus_fourvector.Px()/200)+0.5,(piminus_fourvector.Py()/200)+0.5,0.007,"|>")
    label.DrawPaveLabel((piminus_fourvector.Px()/200)+0.51,(piminus_fourvector.Py()/200)+0.51,(piminus_fourvector.Px()/200)+0.54,(piminus_fourvector.Py()/200)+0.54,"-211")
    c2.Write()
myfile.Write()
final_state_momentum_display()
myfile.Close()
