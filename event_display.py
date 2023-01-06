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


def event_display(directory, name, description):
    
    gEve = TEveManager.Create() 
    
    track_list = TEveTrackList()
    prop = track_list.GetPropagator()
    prop.SetMagField(0)
    track_list.SetElementName(track_list.GetElementName()+ ", zeroB")
    #track = make_track(prop, 1)
    rc = TEveRecTrackD()
    rc.fV.Set(0., 0., 0.)
    rc.fP.Set(1, 1, 1)
    rc.fSign = 1
    track = TEveTrack(rc, prop)
    track.SetName("Charge 1")


    #daughter 0
    pm1 = TEvePathMarkD(TEvePathMarkD.kDaughter)
    pm1.fV.Set(-40., -40., 3.)
    track.AddPathMark(pm1)
    # daughter 1
    pm2 = TEvePathMarkD(TEvePathMarkD.kDaughter)
    pm2.fV.Set(57., -89., -9.)
    track.AddPathMark(pm2)

    track_list.SetLineColor(kMagenta)
    track.SetLineColor(track_list.GetLineColor())
    
    gEve.AddElement(track_list)
    track_list.AddElement(track)
 
    track.MakeTrack()
 
    ev = gEve.GetDefaultViewer()
    gv = ev.GetGLViewer()
    gv.SetGuideState(TGLUtil.kAxesOrigin, kTRUE, kFALSE, 0)
 
    gEve.Redraw3D(kTRUE)
    ROOT.gSystem.ProcessEvents()
 
    #gv.CurrentCamera().RotateRad(-0.5, 1.4)
    gv.RequestDraw()
        
    #gv.SavePicture("plots/event_display0.png")
event_display("data","Initial study","H tau tau events")
myfile.Close()
