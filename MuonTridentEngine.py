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
from Tools import Tools
import BaseCut as BaseCut



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
            print("Alignment pars recovered ")
            self.geo.modules['Scifi'].InitEvent(self.eventTree.EventHeader)
            self.geo.modules['MuFilter'].InitEvent(self.eventTree.EventHeader)
        # if faireventheader, rely on user to select correct geofile.
        self.scifi_hit_density_cut = BaseCut.ScifiHitDensityCut(self.eventTree,self.geo)
        self.scifi_fiducial_cut = BaseCut.FiducialCut()
        self.tools = Tools(self.trackTask)        
 
    def event_delta_t_cut(self,delta_e,delta_t):
        header=self.eventTree.EventHeader
        current_entry = self.eventTree.GetReadEvent()
        current_time  = header.GetEventTime()
        passes = True
        rc= self.eventTree.GetEvent(current_entry+delta_e)
        sign = (delta_e>0)-(delta_e<0)
        if -sign*(current_time-header.GetEventTime()) <= delta_t:
            passes = False
        rc=self.eventTree.GetEvent(current_entry)
        return passes
   

    def filter_time(self):
        filteredHits = ROOT.TClonesArray("sndScifiHit",10000)
        mufilter_hits = ROOT.TClonesArray("MuFilterHit",10000)
        header = self.ioman.GetObject("EventHeader.")    
        filteredHits.BypassStreamer(ROOT.kTRUE)
        self.ioman.Register("EventHeader","sndEventHeader",header,ROOT.kTRUE)
        self.ioman.Register("Digi_ScifiHits", "DigiScifiHit_det", filteredHits, ROOT.kTRUE)
        self.ioman.Register("Digi_MuFilterHits","DigiMuFitlerHit_det",mufilter_hits,ROOT.kTRUE)
        print("things registered")
        B = ROOT.TList()
        branches_to_copy = ["Digi_ScifiHits", "EventHeader","Digi_MuFilterHits"]
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
            if not self.event_delta_t_cut(-1,100):
                continue
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
        flag = ROOT.Event_Type('Event_Type') #ROOT.TClonesArray("TObjString",2) # ROOT.TMap() # ROOT.TClonesArray("TObjString",2) # remember 0 for "scifi" 1 for "ds"
        header = self.ioman.GetObject("EventHeader")
        self.ioman.Register("EventHeader","sndEventHeader",header,ROOT.kTRUE)
        self.ioman.Register("Digi_ScifiHits", "DigiScifiHit_det", scifi_hits, ROOT.kTRUE)
        self.ioman.Register("Digi_MuFilterHits","DigiMuFitlerHit_det",mufilter_hits,ROOT.kTRUE)
        self.ioman.Register("Event_Type","MultiMuonType",flag,ROOT.kTRUE)
        print("things registered")
        B_list = ROOT.TList()
        branches_to_copy = {"Digi_ScifiHits", "EventHeader","Digi_MuFilterHits","Event_Type"}
        B_list.SetName('BranchList')
        for aBranch in branches_to_copy:
            B_list.Add(ROOT.TObjString(aBranch))
        print("branch list was set")
        self.ioman.SetBranchNameList(B_list)
        self.ioman.WriteFolder()
        print("loop started !")
        A,B=ROOT.TVector3(),ROOT.TVector3()
        N_collected = 0
        for n in range(self.eventTree.GetEntries()) :
            self.eventTree.GetEvent(n)
            n3D = {"scifi":[0,0],"ds":[0,0]}
            reco_2dTracks = {"scifi":{0:[],1:[]}, "ds" :{0:[],1:[]}}
            flag.Clear()
#            print(f"Flag Cleared - scifi: {flag.GetValue('scifi')}, ds: {flag.GetValue('ds')}")
            if n%100000==0 : print("Event at ",n)
            rc = self.trackTask.multipleTrackCandidates(nMaxCl=8,dGap=0.2,dMax=0.8,dMax3=0.8,ovMax=1,doublet=True,debug=False)
            for p in range(2):
                for trackId in self.trackTask.multipleTrackStore['trackCand'][p]:
                    if trackId < 1000: continue
                    if trackId in self.trackTask.multipleTrackStore['cloneCand'][p]: continue
                    n3D["scifi"][p]+=1
                    reco_2dTracks["scifi"][p].append(self.trackTask.multipleTrackStore['trackCand'][p][trackId])
            """
            rc=self.trackTask.multipleTrackStore.clear()
            rc = self.trackTask.multipleTrackCandidates2(self.geo.modules['MuFilter'],self.eventTree.Digi_MuFilterHits) 
            for p in range(2):
                for trackId in self.trackTask.multipleTrackStore['trackCand'][p]:
                    if p == 0 and trackId< 100 : continue
                    if p == 1 and trackId< 1000: continue
                    if trackId in self.trackTask.multipleTrackStore['cloneCand'][p]: continue
                    n3D["ds"][p]+=1
                    reco_2dTracks["ds"][p].append(self.trackTask.multipleTrackStore['trackCand'][p][trackId])
            """

            useless_event = {"scifi": True, "ds": True}

# Loop over the keys in n3D
            for key in n3D:
                if n3D[key] == [2, 2]:
                    flag.Add(key,"dimuon")
                    useless_event[key] = False
                elif n3D[key] == [3, 3]:
                    flag.Add(key,"triplemuon")
                    useless_event[key] = False
                elif n3D[key] == [4, 4]:
                    flag.Add(key,"quadramuon")
                    useless_event[key] = False
                elif n3D[key] == [5, 5]:
                    flag.Add(key, "pentamuon")
                    useless_event[key] = False
                else:
                    useless_event[key] = True
# If all events are useless, skip the rest of the loop
            if all(useless_event.values()):
                continue 
            if n3D['scifi']!=[3,3]:continue
            sparse_event = self.scifi_hit_density_cut.pass_cut(reco_2dTracks)
            scifi_fiducial = self.scifi_fiducial_cut.pass_cut(reco_2dTracks)
 
            passes_criteria_scifi = sparse_event and scifi_fiducial
            passes_criteria_ds = False
            if passes_criteria_scifi:
                old_value = flag.GetValue('scifi')
                if passes_criteria_scifi:
                    flag.Add("scifi", 'good_' + old_value)
            else :
                continue





            """
# If ds is missing and scifi is present but doesn't pass criteria, skip the event and vice versa.
            if flag.GetValue('ds') == '' and flag.GetValue('scifi')!='' and not passes_criteria_scifi:
                continue  # Skip if scifi is bad and there's no ds tracks
            if flag.GetValue('scifi')=='' and flag.GetValue('ds')!='' and not passes_criteria_ds:
                continue

# Modify the scifi flag based on whether the criteria were met
            if flag.GetValue('scifi')!='':
                old_value = flag.GetValue('scifi')
                if passes_criteria_scifi:
                    flag.Add("scifi", 'good_' + old_value)
                else:
                    flag.Add("scifi",'bad_' + old_value)

# Modify the ds flag, always adding a "bad" label in this case
            if flag.GetValue('ds')!='' :
                old_value = flag.GetValue('ds')
                if passes_criteria_ds:
                    flag.Add("ds",'good_'+old_value)
                else:    
                    flag.Add('ds','bad_' + old_value)

            
            
# Print the updated flags for scifi and ds
#            print(n3D['scifi'],n3D['ds'])
#            print(f"Writing Event {n}- scifi: {flag.GetValue('scifi')}, ds: {flag.GetValue('ds')}")
            """

            self.ioman.Fill()
            N_collected+=1

        self.ioman.Write()
#        print(N_collected)

