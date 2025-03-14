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
from scipy.optimize import fsolve


#ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gInterpreter.Declare("""
#include "MuFilterHit.h"
#include "AbsMeasurement.h"
#include "TrackPoint.h"

void fixRoot(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllSignals(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRootT(MuFilterHit& aHit,std::vector<int>& key,std::vector<float>& value, bool mask) {
   std::map<int,float> m = aHit.GetAllTimes(false);
   std::map<int, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}
void fixRoot(MuFilterHit& aHit, std::vector<TString>& key,std::vector<float>& value, bool mask) {
   std::map<TString, float> m = aHit.SumOfSignals();
   std::map<TString, float>::iterator it = m.begin();
   while (it != m.end())
    {
        key.push_back(it->first);
        value.push_back(it->second);
        it++;
    }
}

void fixRoot(std::vector<genfit::TrackPoint*>& points, std::vector<int>& d,std::vector<int>& k, bool mask) {
      for(std::size_t i = 0; i < points.size(); ++i) {
        genfit::AbsMeasurement*  m = points[i]->getRawMeasurement();
        d.push_back( m->getDetId() );
        k.push_back( int(m->getHitId()/1000) );
    }
}
""")





class MuonTridentAnalysis:
    "generate histograms for Muon Filter"
    def __init__(self,options):
        print("Initializing Multi Muon Search Engine ... ")
        self.Tkey  = ROOT.std.vector('TString')()
        self.Ikey   = ROOT.std.vector('int')()
        self.Value = ROOT.std.vector('float')()
        self.options=options
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
        if self.f.Get('rawConv'):
            self.eventTree = self.f.rawConv
        else : self.eventTree = self.f.cbmsim
        self.eventTree.GetEvent(0)
        self.run      = ROOT.FairRunAna()
        self.ioman = ROOT.FairRootManager.Instance()
        self.ioman.SetTreeName(self.eventTree.GetName())
        self.outFile = ROOT.TMemFile('dummy','CREATE')
        self.source = ROOT.FairFileSource(self.f)
        self.run.SetSource(self.source)
        self.sink = ROOT.FairRootFileSink(self.outFile)
        if self.f.Get('rawConv') : self.sink.SetRunId(int(options.runNumber))
        else:  self.sink.SetRunId(0)
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
        if self.f.Get("rawConv"): self.run.SetRunId(int(options.runNumber))
        else:  self.run.SetRunId(0)
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
        self.hist={}
        self.systemAndPlanes   = {1:2,2:5}
        self.systemAndBars     = {1:7,2:10}
        self.systemAndChannels = {2:16}
        self.sdict             = {1:'Veto',2:'US'}
        self.cases = {1:"single",2:"dimuon",3:"threemuon"}
        for system in self.systemAndPlanes:
            for plane in range(self.systemAndPlanes[system]):
                if system == 2:
                    ut.bookHist(self.hist,"pass-Exy-"+str(plane)," ",140,-46,-30.,120,35,52)
                    ut.bookHist(self.hist,"total-Exy-"+str(plane)," ",140,-46,-30.,120,35,52)
                    ut.bookHist(self.hist,"residuals-yEx"+str(plane),"Residuals plane"+str(plane), 100, 35,52, 100,-30.,30.)
                    ut.bookHist(self.hist,"residuals-xEx"+str(plane),"Residuals plane"+str(plane), 100, -46,-30, 100,-30.,30.)
                    self.hist["residuals-yEx"+str(plane)].GetXaxis().SetTitle("yEx [cm]")
                    self.hist["residuals-yEx"+str(plane)].GetYaxis().SetTitle("res_y [cm]")
                    self.hist["residuals-xEx"+str(plane)].GetXaxis().SetTitle("xEx [cm]")
                    self.hist["residuals-xEx"+str(plane)].GetYaxis().SetTitle("res_y [cm]")
                    ut.bookHist(self.hist,"nbrHits-"+str(plane)," nbr of hits ",10,0.,10.)
                    self.hist["nbrHits-"+str(plane)].GetXaxis().SetTitle("nbr of hits")
                orientation = self.systemAndOrientation(system,plane)
                for bar in range(self.systemAndBars[system]):
                    bar_key = self.sdict[system]+"plane"+str(plane)+orientation+"bar"+str(bar)
                    if system == 2 or system == 1:
                        for key in range(1,4):
                            ut.bookHist(self.hist, "qdc_"+self.cases[key]+bar_key," ",40,-10,800) # look at the qdc scale !
                for s in range(1,5):            
                    for o in range(2):
                        so = s*10+o
                        if so == 40 : continue
                        ut.bookHist(self.hist,'n_hitsDS_so'+str(so),"",10,0,10)
        self.gr = {}
        dmin_bins = [0,1.,2.5,4.5,7.5,11.5,16.5,25.]  # Variable bin widths in Y
        # Number of bins = (number of edges - 1)
        n_dmin_bins = len(dmin_bins) - 1
        dmin_bins_array = ROOT.std.vector('double')(dmin_bins)
        self.dmax_vs_dmin={}
 

        for p in range(2):
            self.dmax_vs_dmin['data_histo_dmax_vs_dmin_'+str(p)] = ROOT.TH2D("h2"+str(p), "2D Histogram with Variable Y Binning;X axis;Y axis",n_dmin_bins,dmin_bins_array.data(),25,0,25)
            ut.bookHist(self.hist,f'data_dmin_dmax_ratio_{p}',"dmin to dmax ratio",60,0,0.6)
            ut.bookHist(self.hist,"track_dist_"+str(p)," minimum distance between pairs ",100,0,20)
            ut.bookHist(self.hist,"theta_xz","theta xz",100,-0.05,0.05)
            ut.bookHist(self.hist,"theta_yz","theta yz",100,-0.05,0.05)
            ut.bookHist(self.hist,"theta_xz_yz","theta xz vs yz",100,-0.05,0.05,-0.05,0.05)
            ut.bookHist(self.hist,"dmax_vs_dmin_"+str(p),"",50,0,20.,50,0,20.)
            ut.bookHist(self.hist,"single_event_hit_time","",50,0,100)
#            ut.bookCanvas(self.hist,"dmax_vs_dmin"+str(p),"",700,400)
        for p in range(2):
            ut.bookHist(self.hist,"chi2ndof"+str(p),"chi2 ndof",50,0,10)
        
        for s in range(1,5):
            for o in range(2):
                for quality in ["bad","good"]:
                    ut.bookHist(self.hist,'n_hitsDS_so'+str(s*10+o)+quality,'ds hit density' + str(s*10+o),10,0,10)
        ut.bookHist(self.hist,"theta_reco"," ",100,0,200,100,0,0.1)




    def analyze_mu3(self):
        print("loop started !")
        A,B=ROOT.TVector3(),ROOT.TVector3()
        MuFilter= self.geo.modules['MuFilter']
        N_outside = 0
        N_noisy = 0
        dmax = {}
        dmin = {}
        for i in range(2):
            dmax["data_"+str(i)]=np.array([],dtype=float)
            dmin["data_"+str(i)]=np.array([],dtype=float)
        for n in range(self.eventTree.GetEntries()) :
            self.eventTree.GetEvent(n)
            if n%100000==0 : print("Event at ",n)
            rc = self.eventTree.GetEvent(n)
            flag = self.eventTree.Event_Type
            if not flag.GetValue("scifi")=='good_triplemuon' :continue
            n3D = {"scifi":[0,0],"ds":[0,0]}
            reco_2dTracks = {"scifi":{0:[],1:[]}, "ds" :{0:[],1:[]}}
            rc = self.trackTask.multipleTrackCandidates(nMaxCl=8,dGap=0.2,dMax=0.8,dMax3=0.8,ovMax=1,doublet=True,debug=False)
            for p in range(2):
                for trackId in self.trackTask.multipleTrackStore['trackCand'][p]:
                    if trackId < 1000: continue
                    if trackId in self.trackTask.multipleTrackStore['cloneCand'][p]: continue
                    n3D["scifi"][p]+=1
                    reco_2dTracks["scifi"][p].append(self.trackTask.multipleTrackStore['trackCand'][p][trackId])

            rc=self.trackTask.multipleTrackStore.clear()
            rc = self.trackTask.multipleTrackCandidates2(self.geo.modules['MuFilter'],self.eventTree.Digi_MuFilterHits) 
            for p in range(2):
                for trackId in self.trackTask.multipleTrackStore['trackCand'][p]:
                    if p == 0 and trackId< 100 : continue
                    if p == 1 and trackId< 1000: continue
                    if trackId in self.trackTask.multipleTrackStore['cloneCand'][p]: continue
                    n3D["ds"][p]+=1
                    reco_2dTracks["ds"][p].append(self.trackTask.multipleTrackStore['trackCand'][p][trackId])
            acceptance = {'Veto':True,'US':True}
            for p in range(2):
                for aTrack in reco_2dTracks["scifi"][p]:
                    if not self.check_reaching(aTrack,1,0,p) : acceptance['Veto'] = False
                    if not self.check_reaching(aTrack,2,4,p) : acceptance['US'] = False
            for aHit in self.eventTree.Digi_MuFilterHits:
                if aHit.GetSystem() == 3 : continue
                sumSignal = self.map2Dict(aHit,'SumOfSignals')
                detID = aHit.GetDetectorID()
                if not aHit.isValid() : continue
                if not self.is_UShit_isolated(aHit): continue
                system = aHit.GetSystem()
                plane  = (detID%10000)//1000
                bar = detID%1000
                MuFilter.GetPosition(detID,A,B)
                n_passing = 0
                for aTrack in reco_2dTracks["scifi"][0]:
                    rc = aTrack.Fit("pol1","QS")
                    line = aTrack.GetFunction("pol1")
                    ex=line.Eval((A[2]+B[2])/2)
                    if ex<(A[1]+B[1])/2+2.8 and ex>(A[1]+B[1])/2-2.8 : 
                        n_passing+=1
                if n_passing < 1 : continue
                N_fired = len(self.map2Dict(aHit,'GetAllSignals'))
                if N_fired < 6 :
                    if n_passing >1 : N_noisy+=1
                    continue
                if system == 2 and acceptance['US']: 
                    bar_key = "USplane"+str(plane)+"horizontalbar"+str(bar)
                    self.hist["qdc_"+self.cases[n_passing]+bar_key].Fill(sumSignal['Sum'])
                if system == 1 and acceptance['Veto']:
                    bar_key = "Vetoplane"+str(plane)+"horizontalbar"+str(bar)
                    self.hist["qdc_"+self.cases[n_passing]+bar_key].Fill(sumSignal['Sum'])
#            for aHit in self.eventTree.Digi_ScifiHits:
#                print("hit time: " aHit.GetTime())
        


            thetas=[]
            for p in range(2):
                keys = list(reco_2dTracks['scifi'][p])
                pairs={}
                for i in range(len(reco_2dTracks['scifi'][p])-1):
                    for j in range(i+1,len(reco_2dTracks['scifi'][p])):
                        tr_1,tr_2 = keys[i],keys[j]
                        rc = tr_1.Fit("pol1","QS")
                        line1 = tr_1.GetFunction("pol1")
                        rc = tr_2.Fit("pol1","QS")
                        line2 = tr_2.GetFunction("pol1")
                        X1,X2=line1.Eval(280),line2.Eval(280)
                        slope1,slope2=line1.GetParameter(1),line2.GetParameter(1)
                        theta=self.angle_between_lines(slope1, slope2)
                        d=abs(X1-X2)
                        if d not in pairs :
                            pairs[d]=theta
                sorted_dict = dict(sorted(pairs.items()))
                first_key = next(iter(sorted_dict))
                last_key = next(reversed(sorted_dict))
                dmax["data_"+str(p)]=np.append(dmax["data_"+str(p)],last_key)
                dmin["data_"+str(p)]=np.append(dmin["data_"+str(p)],first_key)
                print(last_key,first_key,n)
#                    self.hist["dmax_vs_dmin_"+str(p)].Fill(first_key,last_key)
                dmin_dmax_ratio = first_key/last_key
                self.dmax_vs_dmin['data_histo_dmax_vs_dmin_'+str(p)].Fill(first_key,last_key)
                self.hist[f'data_dmin_dmax_ratio_{p}'].Fill(dmin_dmax_ratio)

            if self.eventTree.EventHeader.isB2noB1():print("bingooo", self.eventTree.GetReadEvent(),self.eventTree.Event_Type.GetValue("scifi"),self.eventTree.Event_Type.GetValue("ds"))
        ut.writeHists(self.hist,"histograms_mu3_qdc"+str(self.options.runNumber)+".root")
        ana_file = ROOT.TFile("histograms_mu3_qdc"+str(self.options.runNumber)+".root","UPDATE")
        ana_file.cd()
        for p in range(2):
            graph_name='data_dmax_vs_dmin_'+str(p)
            self.gr[graph_name]=ROOT.TGraph(len(dmax["data_"+str(p)]),dmin["data_"+str(p)],dmax["data_"+str(p)])        
        for p in range(2):
            graph_name='data_dmax_vs_dmin_'+str(p)
            self.gr[graph_name].Write(graph_name)
            self.dmax_vs_dmin['data_histo_dmax_vs_dmin_'+str(p)].Write('data_histo_dmax_vs_dmin_'+str(p))

        ana_file.Close()
 



    def angle_between_lines(self,slope1, slope2):
        if slope1 == slope2:
            return 0.0  # The lines are parallel
        # Calculate the tangent of the angle between the lines
        tan_theta = abs((slope2 - slope1) / (1 + slope1 * slope2))
    
        # Convert the tangent to an angle in degrees
        angle = ROOT.TMath.ATan(tan_theta) #* (180 / ROOT.TMath.Pi())
    
        return angle





    def map2Dict(self,aHit,T='GetAllSignals',mask=True):
        if T=="SumOfSignals":
            key = self.Tkey
        elif T=="GetAllSignals" or T=="GetAllTimes":
            key = self.Ikey
        else: 
           print('use case not known',T)
           1/0
        key.clear()
        self.Value.clear()
        if T=="GetAllTimes": ROOT.fixRootT(aHit,key,self.Value,mask)
        else:                         ROOT.fixRoot(aHit,key,self.Value,mask)
        theDict = {}
        for k in range(key.size()):
            if T=="SumOfSignals": theDict[key[k].Data()] = self.Value[k]
            else: theDict[key[k]] = self.Value[k]
        return theDict

    def systemAndOrientation(self,s,plane):
        if s==1 or s==2: return "horizontal"
        if plane%2==1 or plane == 6: return "vertical"
        return "horizontal"

    def smallSiPMchannel(self,i):
        if i==2 or i==5 or i==10 or i==13: return True
        else: return False

    def check_reaching(self,theTrack,system,plane,track_projection):
        reaching  = True
        orientation = self.systemAndOrientation(system,plane)
        bars = range(self.systemAndBars[system])
        MuFilter= self.geo.modules['MuFilter']
        A,B=ROOT.TVector3(),ROOT.TVector3()
        if orientation == 'horizontal':
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[0]),A,B)
            x_start, x_end = A[0], B[0]
            y_start = (A[1]+B[1])/2
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[-1]),A,B)
            y_end = (A[1]+B[1])/2
        if orientation == 'vertical':
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[0]),A,B)
            y_start, y_end = A[1], B[1]
            x_start = (A[0]+B[0])/2
            MuFilter.GetPosition(int(system*10**4+plane*1E3+bars[-1]),A,B)
            x_end = (A[0]+B[0])/2
        x_det = [x_start,x_end]
        y_det = [y_start,y_end]
        z_mid = (A[2]+B[2])/2
        ex = theTrack.Eval(z_mid)
        if track_projection == 1 :
            reaching = (min(x_det) < ex and ex < max(x_det))
        if track_projection == 0 :
            reaching = (min(y_det) < ex and ex < max(y_det))
        return reaching

    def is_UShit_isolated(self,theHit):
        detID = theHit.GetDetectorID()
        the_plane  = (detID%10000)//1000
        the_bar = detID%1000
        isolated = True
        the_system = theHit.GetSystem()
        for aHit in self.eventTree.Digi_MuFilterHits:
            if aHit.GetSystem() != the_system : continue
            if aHit.GetDetectorID() == detID : continue 
            if (aHit.GetDetectorID()%10000)//1000 == the_plane:
                if abs(the_bar-aHit.GetDetectorID()%1000)==1:isolated = False
        return isolated
