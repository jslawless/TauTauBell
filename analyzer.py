import pyhepmc as hep
import math,random
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F,TH2F, TFile
from ROOT import gROOT, gBenchmark, TLorentzVector
import os

def muon_boost_vector(muon, antimuon, photons):
    momentum = muon.momentum + antimuon.momentum
    #Doing this to make sure the photons didn't come from
    #the neutral pions, as that would mess up the boosting
    for ph in photons:
        if len(ph.parents) == 0:
            momentum = momentum + ph.momentum
    muon_boost = TLorentzVector()
    muon_boost.SetPtEtaPhiM(momentum.pt(),momentum.eta(),momentum.phi(),125)
    return muon_boost

def spacetime_separation(v1, v2):
    return (v1.t-v2.t)**2 - (v1.x-v2.x)**2 - (v1.y-v2.y)**2 - (v1.z-v2.z)**2

def mock_separation(v1, v2):
    return abs(math.sqrt((v1.x-v2.x)**2 + (v1.y-v2.y)**2 + (v1.z-v2.z)**2)/(v1.t-v2.t))


def polarimeter(tau,pi,nu,neutralpi,no_neutralpion):
    if no_neutralpion:
        return pi.Vect()
    mt = 1.777
    return mt*(pi.Energy() - neutralpi.Energy())*(pi.Vect() - neutralpi.Vect()) + 0.5*((pi + neutralpi).M2())*nu.Vect()
    
    

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
    hist8 = TH1F( name + "_photonEnergy", "Photon Energy", 200,0, 1)
    ssep = TH1F( name + "_spacetimeSeparation", "Spacetime Separation",300, -200, 1)
    ssep.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    fine_ssep = TH1F( name + "_finespacetimeSeparation", "Spacetime Separation",30, -0.5, 1)
    fine_ssep.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    speed = TH1F( name + "_speed", "Speed",300, 0, 20)
    speed.GetXaxis().SetTitle("#Delta r / #Delta t ")
    angle = TH1F( name + "_angle", "Angle Between Polarimeter Vectors, No Neutral Pions", 100, -3.15,3.15)
    angle.GetXaxis().SetTitle("#Delta#phi")
    npangle = TH1F( name + "_npangle", "Angle Between Polarimeter Vectors, Neutral Pions", 100, -3.15,3.15)
    npangle.GetXaxis().SetTitle("#Delta#phi")
    bell_effect = TH2F( name+ "_bellInequality","Spacetime Separation vs Delta Phi Between Unit Vectors", 25, -10,10.0,
                       20, 4, -4)
    bell_effect.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    bell_effect.GetXaxis().SetLimits(0.0,10.0)
    bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")
    mock_bell_effect = TH2F( name+ "_mockBellInequality","Spacetime Separation vs Delta Phi Between Unit Vectors", 25,
                            -10,10,
                       20, 4, -4)
    mock_bell_effect.GetXaxis().SetTitle("#Delta t^{2}-#Delta r^{2} [mm^{2}]")
    mock_bell_effect.GetXaxis().SetLimits(0,10)

    mock_bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")
    
    speed_bell_effect = TH2F( name+ "_speedBellInequality","Speed of Mediator vs Delta Phi Between Unit Vectors", 25, 1,20,20, 4, -4)
    speed_bell_effect.GetXaxis().SetTitle("#Delta r/#Delta t")
    speed_bell_effect.GetXaxis().SetLimits(0.0,10.0)
    speed_bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")
    speed_mock_bell_effect = TH2F( name+ "_speedMockInequality","Speed of Mediator vs Delta Phi Between Unit Vectors", 25,1,20,20, 4, -4)
    speed_mock_bell_effect.GetXaxis().SetTitle("#Delta r/#Delta t")
    speed_mock_bell_effect.GetXaxis().SetLimits(0,10)

    speed_mock_bell_effect.GetYaxis().SetTitle("Angle Between Unit Vectors of the Decay Plane")



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
    count = 0 
    while(not endCondition):
        count = count + 1
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
        tnu_candidates =[]
        antitau_candidates = []
        antitnu_candidates = []
        piplus_candidates = []
        piminus_candidates = []
        photons = []
        pi0s = []
        for candidate in evt.particles:
            if candidate.pid == 13:
                muon_candidates.append(candidate)
            if candidate.pid == -13:
                antimuon_candidates.append(candidate)
            if candidate.pid == 15:
                tau_candidates.append(candidate)
            if candidate.pid == -15:
                antitau_candidates.append(candidate)
            if candidate.pid == 16:
                tnu_candidates.append(candidate)
            if candidate.pid == -16:
                antitnu_candidates.append(candidate)
            if candidate.pid == 211:
                piplus_candidates.append(candidate)
            if candidate.pid == -211:
                piminus_candidates.append(candidate)
            if candidate.pid == 22:
                photons.append(candidate)
            if candidate.pid == 111 or candidate.pid == -111:
                pi0s.append(candidate)
        
        #Find the neutral pions that belong to each tau
        tau_pions = []
        antitau_pions = []
        for i in pi0s:
            if len(i.parents) > 0:
                if i.parents[0].pid == 15:
                    tau_pions.append(i)
                elif i.parents[0].pid == -15:
                    antitau_pions.append(i)
        
        #Only accept events where each tau
        #has one neutral pion or less
        if len(tau_pions) > 1:
            continue
        if len(antitau_pions) > 1:
            continue

        photon_total_energy = 0

        for ph in photons:
            hist8.Fill(ph.momentum.e)
            photon_total_energy = photon_total_energy + ph.momentum.e
        total_photons = total_photons + len(photons)

        muon = muon_candidates[0]
        antimuon = antimuon_candidates[0]
        tau = tau_candidates[0]
        nu = tnu_candidates[0]
        antinu = antitnu_candidates[0]
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
        for i in tnu_candidates:
            if i.momentum.pt() > nu.momentum.pt():
                tau = i
        for i in antitnu_candidates:
            if i.momentum.pt() > antinu.momentum.pt():
                antitau = i
        for i in piplus_candidates:
            if i.momentum.pt() > piplus.momentum.pt():
                piplus = i
        for i in piminus_candidates:
            if i.momentum.pt() > piminus.momentum.pt():
                piminus = i
        momentum = (muon.momentum + antimuon.momentum)
        hist1.Fill(momentum.m())
        
        #Here the different boost vectors are created
        muon_boost = muon_boost_vector(muon,antimuon,photons)
        tau_fourvector = TLorentzVector()
        tau_fourvector.SetPxPyPzE(tau.momentum.px,tau.momentum.py,tau.momentum.pz,tau.momentum.e)
        antitau_fourvector = TLorentzVector()
        antitau_fourvector.SetPxPyPzE(antitau.momentum.px,antitau.momentum.py,antitau.momentum.pz,antitau.momentum.e)
        nu_fourvector = TLorentzVector()
        nu_fourvector.SetPxPyPzE(nu.momentum.px,nu.momentum.py,nu.momentum.pz,nu.momentum.e)
        antinu_fourvector = TLorentzVector()
        antinu_fourvector.SetPxPyPzE(antinu.momentum.px,antinu.momentum.py,antinu.momentum.pz,antinu.momentum.e)
        piplus_fourvector = TLorentzVector()
        piplus_fourvector.SetPxPyPzE(piplus.momentum.px,piplus.momentum.py,piplus.momentum.pz,piplus.momentum.e)
        piminus_fourvector = TLorentzVector()
        piminus_fourvector.SetPxPyPzE(piminus.momentum.px,piminus.momentum.py,piminus.momentum.pz,piminus.momentum.e)
        hist2.Fill((tau.momentum+antitau.momentum).m())
        hist3.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())

        #Boost into the HIggs rest frame
        tau_fourvector.Boost(muon_boost.BoostVector())
        antitau_fourvector.Boost(muon_boost.BoostVector())
        nu_fourvector.Boost(muon_boost.BoostVector())
        antinu_fourvector.Boost(muon_boost.BoostVector())
        piplus_fourvector.Boost(muon_boost.BoostVector())
        piminus_fourvector.Boost(muon_boost.BoostVector())
        hist4.Fill((tau_fourvector + antitau_fourvector).Vect().Mag())

        taupi_fourvector = TLorentzVector()
        if len(tau_pions) > 0:
            taupi_fourvector.SetPxPyPzE(tau_pions[0].momentum.px,tau_pions[0].momentum.py,tau_pions[0].momentum.pz,tau_pions[0].momentum.e)
            taupi_fourvector.Boost(muon_boost.BoostVector())
        antitaupi_fourvector = TLorentzVector()
        if len(antitau_pions) > 0:
            antitaupi_fourvector.SetPxPyPzE(antitau_pions[0].momentum.px,antitau_pions[0].momentum.py,antitau_pions[0].momentum.pz,antitau_pions[0].momentum.e)
            antitaupi_fourvector.Boost(muon_boost.BoostVector())

        #Find the unit vector of the decay plane of the taus
        tau_polarimeter = polarimeter(tau_fourvector,
                                      piminus_fourvector,
                                      nu_fourvector,
                                      taupi_fourvector,
                                      (len(tau_pions)==0)
                                      )
        antitau_polarimeter = polarimeter(antitau_fourvector,
                                          piplus_fourvector,
                                          antinu_fourvector,
                                          antitaupi_fourvector,
                                          (len(antitau_pions)==0)
                                          )

        #The ROOT method .Angle() doesn't have a sign, but we need a two sided
        #distribution to look more sinusoidal. We preserve the sign here and add it back in
        #after .Angle() removes it.
        angleforSign = tau_polarimeter.Cross(antitau_polarimeter).Unit()
        sign = -1
        if (angleforSign.Angle(tau_fourvector.Vect()) < (math.pi/2) ):
            sign = 1
        angle_between_polarimeters = sign*tau_fourvector.Vect().Cross(tau_polarimeter).Unit().Angle(tau_fourvector.Vect().Cross(antitau_polarimeter).Unit())
        if ((len(tau_pions) == 0) and (len(antitau_pions) == 0)):
            angle.Fill(angle_between_polarimeters)        
        else:
            npangle.Fill(angle_between_polarimeters)
        
        sep = spacetime_separation(tau.end_vertex.position,antitau.end_vertex.position)
        mock_sep = mock_separation(tau.end_vertex.position,antitau.end_vertex.position)

        ssep.Fill(sep)
        fine_ssep.Fill(sep)
        speed.Fill(mock_sep)

        rand = random.uniform(-math.pi,math.pi)
        if sep > -25:
            bell_effect.Fill(sep,angle_between_polarimeters)
            if mock_sep < 2:
                mock_bell_effect.Fill(sep,angle_between_polarimeters)
            else:
                mock_bell_effect.Fill(sep,rand)
        if mock_sep < 5:
            speed_bell_effect.Fill(mock_sep,angle_between_polarimeters)
            if mock_sep < 2:
                speed_mock_bell_effect.Fill(mock_sep,angle_between_polarimeters)
            else:
                speed_mock_bell_effect.Fill(mock_sep,rand)
            
    print(count)
    myfile.Write()
    myfile.Close()

histogram("data","cp_0","H tau tau events","cp_phase_0")
histogram("mixing_angle_data","cp_pi_half","H tau tau events","cp_phase_pi_half")
