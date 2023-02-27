import pyhepmc as hep
import math,random
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F,TH2F, TFile
from ROOT import gROOT, gBenchmark, TLorentzVector
import os



def muon_boost_vector(muon, antimuon, photons):
    momentum = muon.momentum + antimuon.momentum
    for ph in photons:
        if len(ph.parents) == 0:
            momentum = momentum + ph.momentum
    muon_boost = TLorentzVector()
    muon_boost.SetPtEtaPhiM(momentum.pt(),momentum.eta(),momentum.phi(),125)
    return muon_boost

def spacetime_separation(v1, v2):
    return (v1.t-v2.t)**2 - (v1.x-v2.x)**2 - (v1.y-v2.y)**2 - (v1.z-v2.z)**2

def histogram(directory, name, description,root_file):
    #Set to UPDATE to leave previous files there
    #Open the .root file in RECREATE mode in order to wipe all histograms previously there
    myfile = TFile( 'RootFiles/'+str(root_file)+'.root', 'RECREATE' )

    hist1 = TH1F( name + "_im", "Invariant Mass of Muon Pairs", 50, 0, 100 )
    hist1.GetXaxis().SetTitle("Mass [GeV]")
    hist2 = TH1F( name + "_tau_im", "Invariant Mass of tau Pairs", 100, 0, 200 )
    hist2.GetXaxis().SetTitle("Mass [GeV]")
    hist3 = TH1F( name + "_combined_fv", "Added Magnitude of Tau Three Vectors", 100, 0, 200 )
    hist3.GetXaxis().SetTitle("Momentum [GeV]")
    hist4 = TH1F( name + "_combined_boosted_fv", "Added Magnitude of Tau Three Vectors Boosted using Muon Frame", 100, 0, 200 )
    hist4.GetXaxis().SetTitle("Momentum [GeV]")
    hist5 = TH1F( name + "_angleBetweenUnitVecsPiPlus", "Angle between Unit vectors of plane, tau cross pi plus", 200,
                 -3.5, 3.5)
    hist5.GetXaxis().SetTitle("#Delta#phi")
    hist6 = TH1F( name + "_angleBetweenTaus", "Angle between Tau Three Vectors", 50, -5,5)
    hist6.GetXaxis().SetTitle("#Delta#phi")
    hist7 = TH1F( name + "_angleBetweenUnitVecsPiMinus", "Angle between Unit vectors of plane, tau cross pi minus", 200,
                 -3.5, 3.5)
    hist7.GetXaxis().SetTitle("#Delta#phi")
    hist8 = TH1F( name + "_photonEnergy", "Photon Energy", 200,0, 1)
    hist9 = TH1F( name + "_spacetimeSeparation", "Spacetime Separation",300, -7, 1)
    hist9.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    hist2D = TH2F( name + "_im_dr", "Invariant Mass of Z vs Tau - AntiTau, no cut", 100,81, 101, 100, 0, 12 )
    hist2D.GetXaxis().SetTitle("Invariant Mass (GeV)")
    hist2D.GetYaxis().SetTitle("Tau - Antitau (GeV)")
    crossproduct2D = TH2F ( name + "_angleBetweenUnitVecs2D","Correlation between the two unit vector angles",
                           200,-4,4,200,-4,4)
    crossproduct2D.GetXaxis().SetTitle("Tau cross Pi Plus")
    crossproduct2D.GetYaxis().SetTitle("Tau cross Pi minus")
    bell_effect = TH2F( name+ "_bellInequality","Spacetime Separation vs Delta Phi Between Unit Vectors", 12, -10,10.0,
                       20, 4, -4)
    bell_effect.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    bell_effect.GetXaxis().SetLimits(0.0,10.0)
    bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")
    mock_bell_effect = TH2F( name+ "_mockBellInequality","Spacetime Separation vs Delta Phi Between Unit Vectors", 12,
                            -10,10,
                       20, 4, -4)
    mock_bell_effect.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    mock_bell_effect.GetXaxis().SetLimits(0,10)

    mock_bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")
    piplus_photon_energy = TH2F( name + "_piplus_photon_energy", "Total Photon Energy vs Angle between Unit vectors of plane, tau cross pi plus", 100,-4, 4, 100, 0, 90 )
    piplus_photon_energy.GetXaxis().SetTitle("Tau cross Pi plus")
    piplus_photon_energy.GetYaxis().SetTitle("Summed Energy of Photons(GeV)")
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
    
    low_energy_photons = 0
    total_photons = 0
    
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
        photons = []
        has_pi0 = False
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
            if candidate.pid == 22:
                photons.append(candidate)
            if candidate.pid == 111 or candidate.pid == -111:
                has_pi0 = True
        if has_pi0:
            continue
        photon_total_energy = 0
        for ph in photons:
            hist8.Fill(ph.momentum.e)
            photon_total_energy = photon_total_energy + ph.momentum.e
        total_photons = total_photons + len(photons)
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
        muon_boost = muon_boost_vector(muon,antimuon,photons)
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

        tau_cross_piplus = tau_fourvector.Vect().Cross(piplus_fourvector.Vect()).Unit()
        antitau_cross_piminus = antitau_fourvector.Vect().Cross(piminus_fourvector.Vect()).Unit()
        tau_cross_piminus = tau_fourvector.Vect().Cross(piminus_fourvector.Vect()).Unit()
        antitau_cross_piplus = antitau_fourvector.Vect().Cross(piplus_fourvector.Vect()).Unit()
        angleforSign1 = tau_cross_piplus.Cross(antitau_cross_piminus).Unit()
        sign1 = -1
        if (angleforSign1.Angle(tau_fourvector.Vect()) < (math.pi/2) ):
            sign1 = 1
        angleforSign2 = tau_cross_piminus.Cross(antitau_cross_piplus).Unit()
        sign2 = -1
        if (angleforSign2.Angle(tau_fourvector.Vect()) < (math.pi/2) ):
            sign2 = 1
        hist5.Fill(sign1*tau_cross_piplus.Angle(antitau_cross_piminus))
        hist6.Fill(tau_fourvector.Vect().Angle(antitau_fourvector.Vect()))
        decay_plane = sign2*tau_cross_piminus.Angle(antitau_cross_piplus)
        hist7.Fill(decay_plane)
        sep = spacetime_separation(tau.end_vertex.position,antitau.end_vertex.position)
        hist9.Fill(sep)
        hist2D.Fill(momentum.m(),(tau_fourvector + antitau_fourvector).Vect().Mag())
        crossproduct2D.Fill(sign2*tau_cross_piminus.Angle(antitau_cross_piplus),sign1*tau_cross_piplus.Angle(antitau_cross_piminus))
        if sep > -10:
            bell_effect.Fill(sep,decay_plane)
        if sep < -10:
           a = 0 
        elif sep < -4.5:
            mock_bell_effect.Fill(sep,random.uniform(-math.pi,math.pi))
        else:
            mock_bell_effect.Fill(sep,decay_plane)
        piplus_photon_energy.Fill(sign1*tau_cross_piplus.Angle(antitau_cross_piminus),photon_total_energy)

    myfile.Write()
    myfile.Close()

histogram("data","cp_0","H tau tau events","cp_phase_0")
histogram("mixing_angle_data","cp_pi_half","H tau tau events","cp_phase_pi_half")
