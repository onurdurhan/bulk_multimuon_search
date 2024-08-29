import ROOT
import rootUtils as ut
import os
import SndlhcMuonReco
import SndlhcGeo
import operator
from rootpyPickler import Unpickler
from array import array
from ROOT import TTree
import numpy as np
from argparse import ArgumentParser
import sys
import gc

ROOT.gInterpreter.ProcessLine('#include "/afs/cern.ch/work/o/onur/SNDLHCSOFT_2024/sw/rhel9_x86-64/sndsw/master-local1/analysis/tools/sndSciFiTools.h"')
#ROOT.gROOT.SetBatch(ROOT.kTRUE)

class MuonTridentEngine:
    "generate histograms for Muon Filter"
    def __init__(self,options):
        print("Initializing Multi Muon Search Engine ... ")
        self.geo = SndlhcGeo.GeoInterface(options.geoFile)
#        self.lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
#        self.lsOfGlobals.Add(self.geo.modules['Scifi'])
#        self.lsOfGlobals.Add(self.geo.modules['MuFilter'])
        #Initialize FairLogger: set severity and verbosity
        self.logger = ROOT.FairLogger.GetLogger()
        self.logger.SetColoredLog(True)
        self.logger.SetLogVerbosityLevel('low')
        self.logger.SetLogScreenLevel('WARNING')
        self.logger.SetLogToScreen(True)
        self.f=ROOT.TFile.Open(options.inputFile)
        self.eventTree = self.f.rawConv
        self.eventTree.GetEvent(0)
        self.run      = ROOT.FairRunAna()
        self.ioman = ROOT.FairRootManager.Instance()
        self.ioman.SetTreeName(self.eventTree.GetName())
#        self.outFile = ROOT.TMemFile('dummy','CREATE')
#        self.outFile="filtered_corr_time.root"
        self.outFile = options.outFile
        self.source = ROOT.FairFileSource(self.f)
        self.run.SetSource(self.source)
        self.sink = ROOT.FairRootFileSink(self.outFile)
        self.sink.SetRunId(int(options.runNumber))
        self.run.SetSink(self.sink)
        self.HT_tasks = {'muon_reco_task_Sf':SndlhcMuonReco.MuonReco(),
                'muon_reco_task_DS':SndlhcMuonReco.MuonReco(),
                'muon_reco_task_nuInt':SndlhcMuonReco.MuonReco()}
        for ht_task in self.HT_tasks.values():
            self.run.AddTask(ht_task)

        import SndlhcTracking
        self.trackTask = SndlhcTracking.Tracking() 
        self.trackTask.SetName('simpleTracking')
        self.run.AddTask(self.trackTask)

        #avoiding some error messages
        self.xrdb = ROOT.FairRuntimeDb.instance()
        self.xrdb.getContainer("FairBaseParSet").setStatic()
        self.xrdb.getContainer("FairGeoParSet").setStatic()

        for ht_task in self.HT_tasks.values():
            #    ht_task.SetParFile(os.environ['SNDSW_ROOT']+"/python/TrackingParams.xml")
            ht_task.SetParFile("/afs/cern.ch/work/o/onur/scripts/MultiMuons/TrackingParams_dimuons.xml")
            ht_task.SetHoughSpaceFormat('linearSlopeIntercept')
            #    ht_task.SetHoughSpaceFormat('normal')
            #    ht_task.SetHoughSpaceFormat('linearIntercepts')
            # force the output of reco task to genfit::Track
            # as the display code looks for such output
            ht_task.ForceGenfitTrackFormat()
        self.HT_tasks['muon_reco_task_Sf'].SetTrackingCase('passing_mu_Sf')
        self.HT_tasks['muon_reco_task_DS'].SetTrackingCase('passing_mu_DS')
        self.HT_tasks['muon_reco_task_nuInt'].SetTrackingCase('nu_interaction_products')
        self.run.SetRunId(int(options.runNumber))
        self.run.Init()
        self.OT = self.sink.GetOutTree()
        self.eventTree = self.ioman.GetInTree()
        self.Reco_MuonTracks = self.trackTask.fittedTracks
        self.Cluster_Scifi   = self.trackTask.clusScifi
        self.eventTree.GetEvent(0)
        if self.eventTree.EventHeader.ClassName() == 'SNDLHCEventHeader':
            self.geo.modules['Scifi'].InitEvent(self.eventTree.EventHeader)
            self.geo.modules['MuFilter'].InitEvent(self.eventTree.EventHeader)
        # if faireventheader, rely on user to select correct geofile.
   

    def filter_time(self):
        filteredHits = ROOT.TClonesArray("sndScifiHit",10000)
        mufilter_hits = ROOT.TClonesArray("MuFilterHit",10000)
#        header  = ROOT.SNDLHCEventHeader()
#        mufilter_hits = ioman.GetObject("Digi_MuFilterHits")
        header = self.ioman.GetObject("EventHeader")
        filteredHits.BypassStreamer(ROOT.kTRUE)
        self.ioman.Register("EventHeader","sndEventHeader,",header,ROOT.kTRUE)
        self.ioman.Register("Digi_ScifiHits", "DigiScifiHit_det", filteredHits, ROOT.kTRUE)
        self.ioman.Register("Digi_MuFilterHits","DigiMuFitlerHit_det",mufilter_hits,ROOT.kTRUE)
        print("things registered")
        B = ROOT.TList()
        branches_to_copy = {"Digi_ScifiHits", "EventHeader","Digi_MuFilterHits"}
        B.SetName('BranchList')
        for aBranch in branches_to_copy:
            B.Add(ROOT.TObjString(aBranch))
        print("branch list was set")
        self.ioman.SetBranchNameList(B)
        self.ioman.WriteFolder()
        print("loop started !")
        for n in range(self.eventTree.GetEntries()) :
            self.eventTree.GetEvent(n)
            if n%100000==0 : print("Event at ",n)
            filteredHits.Clear("C")
            mufilter_hits.Clear("C") 
            digi_scifi_hits = self.eventTree.Digi_ScifiHits
            if digi_scifi_hits.GetEntries() < 1 : continue
            hits_in_all_planes = True
            for s in range(1,6):
                for o in range(2):
                    if len(ROOT.snd.analysis_tools.getScifiHits(digi_scifi_hits,s,o)) < 1 :
                        hits_in_all_planes = False
            if not hits_in_all_planes : continue
#            print(len( digi_scifi_hits),n,self.eventTree.Digi_MuFilterHits.GetEntries())
            index = 0
            filtered_hit_container =  ROOT.snd.analysis_tools.filterScifiHits(digi_scifi_hits , 0, "TI18")
            for aHit in filtered_hit_container:
#                print(header.GetRunId(),eventTree.EventHeader.GetRunId())
                filteredHits[index]=aHit
                index+=1
            index = 0
#            if len(digi_scifi_hits) > len(filteredHits) : print("hi at ", n , len(digi_scifi_hits) , len(filteredHits))
            for aHit in self.eventTree.Digi_MuFilterHits:
                mufilter_hits[index]=aHit
                index+=1
            self.ioman.Fill()
        self.ioman.Write()

    def search_mu3(self):
        scifi_hits = self.ioman.GetObject('Digi_ScifiHits')
        mufilter_hits = self.ioman.GetObject("Digi_MuFilterHits")
        header = self.ioman.GetObject("EventHeader")
        self.ioman.Register("EventHeader","sndEventHeader,",header,ROOT.kTRUE)
        self.ioman.Register("Digi_ScifiHits", "DigiScifiHit_det", scifi_hits, ROOT.kTRUE)
        self.ioman.Register("Digi_MuFilterHits","DigiMuFitlerHit_det",mufilter_hits,ROOT.kTRUE)
        print("things registered")
        B_list = ROOT.TList()
        branches_to_copy = {"Digi_ScifiHits", "EventHeader","Digi_MuFilterHits"}
        B_list.SetName('BranchList')
        for aBranch in branches_to_copy:
            B_list.Add(ROOT.TObjString(aBranch))
        print("branch list was set")
        self.ioman.SetBranchNameList(B_list)
        self.ioman.WriteFolder()
        print("loop started !")
        A,B=ROOT.TVector3(),ROOT.TVector3()
        for n in range(self.eventTree.GetEntries()) : #      self.eventTree.GetEntries()):
            self.eventTree.GetEvent(n)
            if n%100000==0 : print("Event at ",n)
            busy_event = {0:False,1:False}
            rc = self.eventTree.GetEvent(n)
            rc = self.trackTask.multipleTrackCandidates(nMaxCl=8,dGap=0.2,dMax=0.8,dMax3=0.8,ovMax=1,doublet=True,debug=False)
            sorted_clusters = {s * 10 + o: [] for s in range(1, 6) for o in range(2)}
            for aCl in self.trackTask.clusScifi:
                so = aCl.GetFirst()//100000
                sorted_clusters[so].append(aCl)
            n3D = [0,0]
            bad_reco = False
            for p in range(2):
                for trackId in self.trackTask.multipleTrackStore['trackCand'][p]:
                    if trackId < 100000 and not self.trackTask.multipleTrackStore['doublet']: continue
                    if trackId in self.trackTask.multipleTrackStore['cloneCand'][p]: continue
                    n3D[p]+=1
                    for s in range(1,6):
                        mat  = 1
                        sipm = 1
                        channel = 64
                        plane_id=channel+1000*sipm+10000*mat+100000*p+1000000*s
                        self.geo.modules['Scifi'].GetSiPMPosition(plane_id,A,B)
                        ex = self.trackTask.multipleTrackStore['trackCand'][p][trackId].Eval((A[2]+B[2])/2)
                        rho = 0
                        for aHit in self.eventTree.Digi_ScifiHits:
                            if not aHit.isValid():continue
                            detID = aHit.GetDetectorID()
                            so = detID//100000
                            self.geo.modules['Scifi'].GetSiPMPosition(aHit.GetDetectorID(),A,B)
                            if so != s*10+p: continue
                            if p == 0 : pos=(A[1]+B[1])/2
                            else : pos=(A[0]+B[0])/2
                            if abs(pos-ex)<1 : rho+=1
                        if rho>10:busy_event[p]=True
                    rc = self.trackTask.multipleTrackStore['trackCand'][p][trackId].Fit("pol1","QS")
                    ndof=self.trackTask.multipleTrackStore['trackCand'][p][trackId].GetFunction("pol1").GetNDF()
                    chi2=self.trackTask.multipleTrackStore['trackCand'][p][trackId].GetFunction("pol1").GetChisquare() 
                    if chi2/ndof>0.05 : bad_reco = True
            if n3D != [3,3] : continue
            if self.eventTree.EventHeader.isB2noB1():print("bingooo")
            if busy_event[0] and busy_event[1] : continue 
            if bad_reco :  continue
            self.ioman.Fill()
        self.ioman.Write()
