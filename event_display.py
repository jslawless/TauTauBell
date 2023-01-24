import pyhepmc as hep
import math
import scipy.constants as sci
from ROOT import TCanvas, TPad, TH1F,TH2F, TFile
from ROOT import gROOT, gBenchmark, TLorentzVector
from ROOT import TEveWindow, TEveWindowManager, TEveManager, TEveTrack,TEveTrackPropagator, TEveViewer, TGLViewer, TEveTrackList, TEveViewer, TEveRecTrackD
from ROOT import TGLUtil, TSystem, TEvePathMarkD
from ROOT import kMagenta, kTRUE, kFALSE
import ROOT
import os


#Set to UPDATE to leave previous files there
#Open the .root file in RECREATE mode in order to wipe all histograms previously there
myfile = TFile( 'RootFiles/event_display.root', 'RECREATE' )


def event_display():
   

    gEve = TEveManager.Create() 
    
    track_list = TEveTrackList()
    prop = track_list.GetPropagator()
    prop.SetMagField(0)

    gEve.AddElement(track_list)

    track_list.SetLineColor(kMagenta)
    reader = hep.ReaderAscii('./single_event.hepmc')

    events = []
    evt = hep.GenEvent()
    endCondition = False
    reader.read_event(evt)
    i = 0.

    for particle in evt.particles: 

        start = particle.production_vertex
        end = particle.end_vertex

        rc = TEveRecTrackD()
        rc.fV.Set(i, i, i)
        rc.fP.Set(10.0*(i+1), 10.0*(i+1), 10*(i+1))
        rc.fSign = 1
        track = TEveTrack(rc, prop)
        track.SetName("Particle " + str(i))

        track_list.AddElement(track)

        track.SetLineColor(track_list.GetLineColor())
        track.MakeTrack()
        print(i)
        i = i+1
    ev = gEve.GetDefaultViewer()
    gv = ev.GetGLViewer()
    gv.SetGuideState(TGLUtil.kAxesOrigin, kTRUE, kFALSE, 0)
 
    gEve.Redraw3D(kTRUE)
    ROOT.gSystem.ProcessEvents()
 
    #gv.CurrentCamera().RotateRad(-0.5, 1.4)
    gv.RequestDraw()
        
    #gv.SavePicture("plots/event_display0.png")
event_display()
myfile.Close()
