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
import json


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file data and MC",default="",required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", default=os.environ["EOSSHIP"]+"/eos/experiment/sndlhc/convertedData/physics/2022/geofile_sndlhc_TI18_V0_2022.root")
parser.add_argument("-j","--jobID",dest="jobID",help="job ID",default=0,required=False)
options = parser.parse_args()


#ROOT.gROOT.SetBatch(ROOT.kTRUE)

ROOT.gInterpreter.Declare("""
        #include "Hit2MCPoints.h"
void fixRoot(Hit2MCPoints& point,std::vector<int>& key,std::vector<float>& value, int detID) {
        std::unordered_map m = point.wList(detID);
        std::unordered_map<int, float>::iterator it = m.begin();;
        while (it != m.end())
        {
            key.push_back(it->first);
            value.push_back(it->second);
            it++;
            }
    }
""")

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





#f_name="/eos/experiment/sndlhc/MonteCarlo/MuonBackground/multi_muons_down/scoring_1.8_Bfield/sndLHC.Ntuple-TGeant4_MuMupairProduction_2muons_3e7pr_digCPP.root"

f_name = options.inputFile
print(f_name,"is the file" )


geo_file = "~/geofile_full.Ntuple-TGeant4_boost1000.0.root"

geo = SndlhcGeo.GeoInterface(geo_file)
lsOfGlobals = ROOT.gROOT.GetListOfGlobals()
lsOfGlobals.Add(geo.modules['Scifi'])
lsOfGlobals.Add(geo.modules['MuFilter'])


#Initialize FairLogger: set severity and verbosity
logger = ROOT.FairLogger.GetLogger()
logger.SetColoredLog(True)
logger.SetLogVerbosityLevel('low')
logger.SetLogScreenLevel('WARNING')
logger.SetLogToScreen(True)

run      = ROOT.FairRunAna()
ioman = ROOT.FairRootManager.Instance()

f=ROOT.TFile.Open(f_name)

eventTree = f.cbmsim
runId = 'sim'

outFile = ROOT.TMemFile('dummy','CREATE')
source = ROOT.FairFileSource(f)
run.SetSource(source)
sink = ROOT.FairRootFileSink(outFile)
run.SetSink(sink)

HT_tasks = {'muon_reco_task_Sf':SndlhcMuonReco.MuonReco(),
            'muon_reco_task_DS':SndlhcMuonReco.MuonReco(),
            'muon_reco_task_nuInt':SndlhcMuonReco.MuonReco()}
for ht_task in HT_tasks.values():
    run.AddTask(ht_task)

import SndlhcTracking
trackTask = SndlhcTracking.Tracking() 
trackTask.SetName('simpleTracking')
run.AddTask(trackTask)

#avoiding some error messages
xrdb = ROOT.FairRuntimeDb.instance()
xrdb.getContainer("FairBaseParSet").setStatic()
xrdb.getContainer("FairGeoParSet").setStatic()

for ht_task in HT_tasks.values():
#    ht_task.SetParFile(os.environ['SNDSW_ROOT']+"/python/TrackingParams.xml")
    ht_task.SetParFile("/afs/cern.ch/user/o/onur/TrackingParams_dimuons.xml")
    ht_task.SetHoughSpaceFormat('linearSlopeIntercept')
#    ht_task.SetHoughSpaceFormat('normal')
#    ht_task.SetHoughSpaceFormat('linearIntercepts')


    # force the output of reco task to genfit::Track
    # as the display code looks for such output
    ht_task.ForceGenfitTrackFormat()
HT_tasks['muon_reco_task_Sf'].SetTrackingCase('passing_mu_Sf')
HT_tasks['muon_reco_task_DS'].SetTrackingCase('passing_mu_DS')
HT_tasks['muon_reco_task_nuInt'].SetTrackingCase('nu_interaction_products')

run.Init()
OT = sink.GetOutTree()
eventTree = ioman.GetInTree()
Reco_MuonTracks = trackTask.fittedTracks
Cluster_Scifi         = trackTask.clusScifi
eventTree.GetEvent(0)
if eventTree.EventHeader.ClassName() == 'SNDLHCEventHeader':
   geo.modules['Scifi'].InitEvent(eventTree.EventHeader)
   geo.modules['MuFilter'].InitEvent(eventTree.EventHeader)
# if faireventheader, rely on user to select correct geofile.

if eventTree.GetBranch('Digi_MuFilterHit'):
# backward compatbility for early converted events
  eventTree.GetEvent(0)
  eventTree.Digi_MuFilterHits = eventTree.Digi_MuFilterHit



class MuonTridentEngine:
    "generate histograms for Muon Filter"
    def __init__(self):
        print("Initializing Muon Filter Engine ...")
        self.Tkey  = ROOT.std.vector('TString')()
        self.Ikey   = ROOT.std.vector('int')()
        self.Value = ROOT.std.vector('float')() 
        self.hist={}
        self.trident_type = {'MuonToMuonPair':0,'GammaToMuonPair':0,'PositronToMuonPair':0}
        ut.bookHist(self.hist,"stage1_pass"," ",25,0,2000.)
        ut.bookHist(self.hist,"stage1_all"," ",25,0,2000.)
        ut.bookHist(self.hist,"stage2_pass"," ",10,0,1.)
        ut.bookHist(self.hist,"stage2_total"," ",10,0,1.)
        ut.bookHist(self.hist,"xz","",330,-6000.,600.,200,-200.,200.)
        ut.bookHist(self.hist,"yz","",330,-6000.,600.,200,-200.,200.)
        ut.bookHist(self.hist,"xy","",200,-300,300,200,-300,300)
        ut.bookHist(self.hist,"z_int","",640,-6000,400)
        for o in range(2):
            for s in range(1,6):
                ut.bookHist(self.hist,"n_clu_so"+str(s*10+o)+"_MC","cluster density "+str(s*10+o),40,0,40)
        self.p1, self.p2 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        for key in self.trident_type:
            ut.bookHist(self.hist,"inv_"+key,"",25,0.,1.)
            ut.bookHist(self.hist,"daughters_mom_"+key,"",25,0,2000.)  
            if key == 'GammaToMuonPair':
                ut.bookHist(self.hist,key+'_theta','',100,0.,1.)
            else : ut.bookHist(self.hist,key+'_theta','',100,0.,0.05)
            ut.bookHist(self.hist,"coplanarity_"+key," ",100,0.,0.05)
            ut.bookHist(self.hist,"xz_dist_"+key,"",25,0,25)
            ut.bookHist(self.hist,"yz_dist_"+key,"",25,0,25)
        ut.bookHist(self.hist,"theta_xz","theta xz",100,0,0.4)
        ut.bookHist(self.hist,"theta_yz","theta yz",100,0,0.4)
        ut.bookHist(self.hist,"theta_xz_yz","theta xz vs yz",100,0.,0.5)
        for p in range(2):
            ut.bookHist(self.hist,"chi2ndof"+str(p),"chi2 ndof",100,0,0.3) 
        ut.bookHist(self.hist,"dz_reco"," reco z res", 100,-1000,1000)
        self.systemAndPlanes   = {2:5}
        self.systemAndBars     = {2:10}
        self.systemAndChannels = {2:16}
        self.sdict             = {2:'US'}
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
                    if system == 2:
                        for key in range(1,4):
                            ut.bookHist(self.hist, "qdc_"+self.cases[key]+bar_key," ",30,0.,300)
                            print("qdc_"+self.cases[key]+bar_key)
       

        self.tools = Tools()
        self.cuts = []
        for cut in BaseCut.__subclasses__():
            print("cuts are", cut)
            self.cuts.append(cut())

#    check_target = lambda aVx : all(b[0] < p < b[1] for p, b in zip(aVx, [(-47.7758,-5.5467),(13.0335 ,55.2793),(289.3743,353.3114)]))

    def efficiency(self):
        proc_names = ['MuonToMuonPair','GammaToMuonPair','PositronToMuonPair']
        weighted = {}
        unweighted = {}
        cases = {1:"single",2:"dimuon",3:"threemuon"}
        for proc in proc_names:
            weighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0}
            unweighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0}
        p1, p2 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        check_target = lambda aVertex : all(b[0] < p < b[1] for p, b in zip(aVertex, [(-47.7758,-5.5467),(13.0335 ,55.2793),(-1000000,288.)])) ###whole wall : 353.3114
        print("Total nbr of events in this file, ",eventTree.GetEntries())
        for n in range(eventTree.GetEntries()):
            if n%1000==0: print('event at %d' % n)
            eventTree.GetEvent(n)
            vertices = self.vertex_array(eventTree.MCTrack) 
            w = eventTree.MCTrack[0].GetWeight()
            h2p = eventTree.Digi_ScifiHits2MCPoints.First()
            multiple_entries = 0
            clusters = self.tools.scifi_clusters()
            tracks_counted = False
            if len(vertices)>1 : print("more than one vertices")
            for aVx in vertices:
#                if not check_target(aVx.get_vertex_pos()) : continue
                self.hist["xz"].Fill(aVx.get_vertex_pos()[2],aVx.get_vertex_pos()[0])
                self.hist["yz"].Fill(aVx.get_vertex_pos()[2],aVx.get_vertex_pos()[1])
                self.hist["xy"].Fill(aVx.get_vertex_pos()[0],aVx.get_vertex_pos()[1])
                self.hist["z_int"].Fill(aVx.get_vertex_pos()[2])
#                self.hist["inv_"+aVx.get_trident_type()].Fill()
#                self.hist["daughters_mom_"+aVx.get_triedent_type()].Fill()
#                self.hist[aVx.get_trident_type()+"_theta"].Fill()
                if aVx.get_vertex_pos()[2]>280. : continue 
                daughters = aVx.get_daughters()
                multiple_entries+=1
                unweighted[aVx.get_trident_type()]['norm']+=1
                weighted[aVx.get_trident_type()]['norm']+=w
                if len(daughters)!=2 : continue
                three_clu = self.three_clusters(h2p, daughters[0], daughters[1], clusters)
                if self.check_three_accepted_clusterwise(three_clu):
#                if self.check_three_accepted(h2p,daughters[0],daughters[1]):
                    unweighted[aVx.get_trident_type()]['inAcceptance']+=1
                    weighted[aVx.get_trident_type()]['inAcceptance']+=w
                    if self.with_original_muon(three_clu):
                        unweighted[aVx.get_trident_type()]['with_original_muon']+=1
                        weighted[aVx.get_trident_type()]['with_original_muon']+=w
                    else : continue
                    if self.check_all_planes(clusters): #self.check_planes_with_hits() :#self.check_all_planes(clusters):
                        unweighted[aVx.get_trident_type()]['hits_exist_in_all']+=1
                        weighted[aVx.get_trident_type()]['hits_exist_in_all']+=w 
                        if self.recoble(three_clu):
                            unweighted[aVx.get_trident_type()]['recoble']+=1
                            weighted[aVx.get_trident_type()]['recoble']+=w 
                            if tracks_counted : continue
#                            trackTask.multipleTrackCandidates(nMaxCl=8,dGap=0.2,dMax=0.8,dMax3=0.8,ovMax=1,doublet=True,debug=False)
                            trackTask.multipleTrackCandidates2(geo.modules['MuFilter'],eventTree.Digi_MuFilterHits)
                            n3D = [0,0]
                            for p in range(2):
                                for trackId in trackTask.multipleTrackStore['trackCand'][p]:
                                    if trackId < 10000 and not trackTask.multipleTrackStore['doublet']: continue
                                    if trackId in trackTask.multipleTrackStore['cloneCand'][p]: continue
                                    n3D[p]+=1
############### DONT FORGET ##################3
                            print(n3D, eventTree.GetReadEvent())
                            continue
                            if  n3D == [3,3] :
                                unweighted[aVx.get_trident_type()]['3/3']+=1
                                weighted[aVx.get_trident_type()]['3/3']+=w
                            if  n3D == [2,3] or n3D == [3,2]:
                                unweighted[aVx.get_trident_type()]['2/3']+=1
                                weighted[aVx.get_trident_type()]['2/3']+=w
                            if n3D == [2,2]:
                                unweighted[aVx.get_trident_type()]['2/2']+=1
                                weighted[aVx.get_trident_type()]['2/2']+=w
#                                print("a weird event detected at event ", n, "the nbr of interactions in the event :", len( vertices))
                            tracks = {0:[],1:[]}
                            if n3D == [3,3]:
                                sparse = {0:True, 1:True}
                                for p in range(2):
                                    keys = list(trackTask.multipleTrackStore['trackCand'][p].keys())
                                    for trackId in trackTask.multipleTrackStore['trackCand'][p]:
                                        if trackId < 100000 and not trackTask.multipleTrackStore['doublet']: continue
                                        if trackId in trackTask.multipleTrackStore['cloneCand'][p]: continue
                                        sparse[p]=self.sparse_event(clusters,trackTask.multipleTrackStore['trackCand'][p][trackId],p,w)
                                        rc=trackTask.multipleTrackStore['trackCand'][p][trackId].Fit("pol1","QS")
                                        theta = trackTask.multipleTrackStore['trackCand'][p][trackId].GetFunction("pol1").GetParameter(1)
                                        tracks[p].append(trackId)
                                        if p == 0 :
                                            self.hist["theta_xz"].Fill(theta)
                                        if p == 1:
                                            self.hist["theta_yz"].Fill(theta)
                                        ndof=trackTask.multipleTrackStore['trackCand'][p][trackId].GetFunction("pol1").GetNDF()
                                        chi2=trackTask.multipleTrackStore['trackCand'][p][trackId].GetFunction("pol1").GetChisquare() 
                                        self.hist["chi2ndof"+str(p)].Fill(chi2/ndof,w*40*2.58)
#                                    print(f"The vertex position is: {vertex} at projection {p} true is vertex {true_pos}")
                                if sparse[0] or sparse[1] : 
                                    unweighted[aVx.get_trident_type()]['sparse']+=1
                                    weighted[aVx.get_trident_type()]['sparse']+=w
                                    distances = {0:[],1:[]}
                                    for p in range(2):
                                        keys = list(trackTask.multipleTrackStore['trackCand'][p].keys())
                                        for i in range(len(trackTask.multipleTrackStore['trackCand'][p])-1):
                                            for j in range(i+1,len(trackTask.multipleTrackStore['trackCand'][p])):
                                                tr_1,tr_2 = keys[i],keys[j]
                                                if tr_1 < 100000 and not trackTask.multipleTrackStore['doublet']: continue
                                                if tr_1 in trackTask.multipleTrackStore['cloneCand'][p]: continue
                                                if tr_2 < 100000 and not trackTask.multipleTrackStore['doublet']: continue
                                                if tr_2 in trackTask.multipleTrackStore['cloneCand'][p]: continue
                                                rc=trackTask.multipleTrackStore['trackCand'][p][tr_1].Fit("pol1","QS")
                                                pos1 = trackTask.multipleTrackStore['trackCand'][p][tr_1].GetFunction("pol1").GetParameter(0)
                                                rc = trackTask.multipleTrackStore['trackCand'][p][tr_2].Fit("pol1","QS")
                                                pos2 = trackTask.multipleTrackStore['trackCand'][p][tr_2].GetFunction("pol1").GetParameter(0)
                                                dist=abs(pos1-pos2)
                                                distances[p].append(dist)
                                        if p == 0: self.hist["xz_dist_"+aVx.get_trident_type()].Fill(max(distances[p]),w*18.959999999999997)
                                        if p == 1: self.hist["yz_dist_"+aVx.get_trident_type()].Fill(max(distances[p]),w*18.959999999999997)
 
                            tracks_counted = True
                            if len(vertices)>1 : continue
                            for aHit in eventTree.Digi_MuFilterHits:
                                detID = aHit.GetDetectorID()
                                the_dict=self.map2Dict(eventTree.Digi_MuFilterHits2MCPoints.First(),detID)
#                                print("the_dict = ",the_dict)
                                n_3 = []
                                for key in the_dict:
                                    trackID = eventTree.MuFilterPoint[key].GetTrackID()
                                    if abs(eventTree.MCTrack[trackID].GetPdgCode())==13:
                                        n_3.append(trackID)
                                if len(n_3)<1 or len(n_3)>3:continue
                                system = aHit.GetSystem()
                                if system != 2 :continue
                                sumSignal = self.map2Dict_qdc(aHit,'SumOfSignals')
                                plane  = (detID%10000)//1000
                                bar = detID%1000
                                bar_key = "USplane"+str(plane)+"horizontalbar"+str(bar)
                                self.hist["qdc_"+cases[len(n_3)]+bar_key].Fill(sumSignal['Sum'])


 
                                
        print("UNWEIGHTED ", unweighted)
        print("WEIGHTED: ", weighted)

        # Write dictionary to a file in JSON format
        with open('data_unweighted_'+str(options.jobID)+'.json', 'w') as json_file:
            json.dump(unweighted, json_file)
        with open('data_weighted_'+str(options.jobID)+'.json','w') as json_file:
            json.dump(weighted, json_file)

#        with open('data'+str(options.jobID)+'.json', 'r') as json_file:
#            loaded_dict = json.load(json_file)
#            print(loaded_dict)
        ut.writeHists(self.hist,"histos_newdigi"+str(options.jobID)+".root")
        return True
        if False :    
##################### SKIP THE REST FOR THE MOMENT ###########################
                theta = self.calculate_angle(daughters[0],daughters[1])
                self.hist["stage2_total"].Fill(theta,w)
                for cut in self.cuts:
                    if cut.pass_cut():
                        self.hist["stage1_pass"].Fill(eventTree.MCTrack[0].GetP(),w)
                        self.hist["stage2_pass"].Fill(theta,w)

                for p in eventTree.ScifiPoint:
                    break
                    if p.GetTrackID()==0 and p.GetZ()>aVx.get_vertex_pos()[2]:
                        px_o, py_o, pz_o = p.GetPx(), p.GetPy(), p.GetPz()
                        break
                for p in eventTree.MuFilterPoint:
                    break
                    if p.GetTrackID()==0 and p.GetZ()<aVx.get_vertex_pos()[2]:
                        px_ref, py_ref, pz_ref = p.GetPx(), p.GetPy(), p.GetPz()
                        break
#        ut.writeHists(self.hist, "efficiency-histos"+f_name[-11:])
        return self.hist
    def is_UShit_isolated(self,theHit):
        detID = theHit.GetDetectorID()
        the_plane  = (detID%10000)//1000
        the_bar = detID%1000
        isolated = True
        for aHit in self.eventTree.Digi_MuFilterHits:
            if aHit.GetSystem()!=2 : continue
            if aHit.GetDetectorID() == detID : continue 
            if (aHit.GetDetectorID()%10000)//1000 == the_plane:
                if abs(the_bar-aHit.GetDetectorID()%1000)==1:isolated = False
        return isolated

    def systemAndOrientation(self,s,plane):
        if s==1 or s==2: return "horizontal"
        if plane%2==1 or plane == 6: return "vertical"
        return "horizontal"


    def map2Dict_qdc(self,aHit,T='GetAllSignals',mask=True):
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


    
    def sparse_event(self,clusters,track,p,w):
        sorted_clusters=self.tools.scifi_clusters_per_plane(clusters)
        sparse = True
        A,B = ROOT.TVector3(),ROOT.TVector3()
        for so in sorted_clusters:
            if not so%10 == p :continue
            aCl = sorted_clusters[so][0]
            aCl.GetPosition(A,B)
            ex = track.Eval(A[2])
            rho=0
            for aCl in sorted_clusters[so]:
                aCl.GetPosition(A,B)
                if p==0 : pos=(A[1]+B[1])/2
                else : pos=(A[0]+B[0])/2
                if abs(pos-ex)<1 : rho+=aCl.GetN()
            self.hist["n_clu_so"+str(so)+"_MC"].Fill(rho,w*40*2.58)
            if rho>10: sparse=False
        return sparse
 


    def clusterType(self,aCl,h2p):
        DetID2Key = {}
        for nHit in range(eventTree.Digi_ScifiHits.GetEntries()):
            DetID2Key[eventTree.Digi_ScifiHits[nHit].GetDetectorID()] = nHit
        for nh in range(aCl.GetN()):
            detID = aCl.GetFirst() # + nh
            aHit = eventTree.Digi_ScifiHits[DetID2Key[detID]]
            theDict = self.map2Dict(h2p,detID)
        return theDict

    def three_clusters(self,h2p,daughter1,daughter2,clusters):
        sorted_clusters = self.tools.scifi_clusters_per_plane(clusters)
        three_clu = {k : {s * 10 + o: [] for s in range(1, 6) for o in range(2)} for k in range(3)}
        DetID2Key = {}
        for nHit in range(eventTree.Digi_ScifiHits.GetEntries()):
            DetID2Key[eventTree.Digi_ScifiHits[nHit].GetDetectorID()] = nHit
        for so in sorted_clusters:
            for aCl in sorted_clusters[so]:
                for nh in range(aCl.GetN()):
                    detID = aCl.GetFirst()+nh
                    aHit = eventTree.Digi_ScifiHits[DetID2Key[detID]]
                    theDict = self.map2Dict(h2p,detID)
                    for link_point in theDict:
                        if eventTree.ScifiPoint[link_point].GetTrackID()< -1 :continue
                        registered_track = eventTree.MCTrack[eventTree.ScifiPoint[link_point].GetTrackID()]
                        if registered_track == daughter1 and aCl not in three_clu[1][so] : three_clu[1][so].append(aCl)
                        if registered_track == daughter2 and aCl not in three_clu[2][so] : three_clu[2][so].append(aCl)
                        if registered_track == eventTree.MCTrack[0] and aCl not in three_clu[0][so] : three_clu[0][so].append(aCl)
        return three_clu

    def check_three_accepted_clusterwise(self, three_clu):
        accepted={0:True,1:True,2:True}
        for k in three_clu:
            nCl_hor = 0
            nCl_ver = 0
            for so in three_clu[k]:
                if so%2==1:
                    nCl_ver+=len(three_clu[k][so])
                else:
                    nCl_hor+=len(three_clu[k][so])
            if nCl_ver< 1 or nCl_hor < 1 : accepted[k] = False
        return all(accepted.values())

    def with_original_muon(self,three_clu):
        n_original_mu = 0
        for so in three_clu[0]:
            n_original_mu+=len(three_clu[0][so])
        if n_original_mu<1: return False
        else : return True


    def check_three_accepted(self,h2p, daughter1, daughter2):
        three_hits = {k : False for k in range(3)}
        for aHit in eventTree.Digi_ScifiHits:
            detID = aHit.GetDetectorID()
            theDict = self.map2Dict(h2p,detID)
            for link_point in theDict:
                registered_track = eventTree.MCTrack[eventTree.ScifiPoint[link_point].GetTrackID()]
                if registered_track == daughter1 : three_hits[1] = True
                if registered_track == daughter2 : three_hits[2] = True
                if registered_track == eventTree.MCTrack[0]: three_hits[0] = True
                
        return all(three_hits.values())

 

    def check_planes_with_hits(self):
        planes_with_hits = {s * 10 + o: False for s in range(1, 6) for o in range(2)}
        for aHit in eventTree.Digi_ScifiHits :
            detID = aHit.GetDetectorID()
            s = int(detID/1000000)
            o = int(detID/100000)%10
            so = 10*s+o
            planes_with_hits[so]=True
        return all(planes_with_hits.values())

    def check_all_planes(self,clusters):
        sorted_clusters = self.tools.scifi_clusters_per_plane(clusters)
        for so in sorted_clusters:
            if len(sorted_clusters[so]) < 1 :return False
        return True

    def recoble(self,three_clu):
        planes = [s * 10 + o for s in range(1, 6) for o in range(2)]
        A,B = ROOT.TVector3(),ROOT.TVector3()
        nOverlap = {'12':0,'13':0,'23':0}
        for k in three_clu:
            planes_with_hit_hor = 0
            planes_with_hit_ver = 0
            for so in three_clu[k]:
                if so%2==1:
                    if len(three_clu[k][so])>0 : planes_with_hit_ver += 1
                else:
                    if len(three_clu[k][so])>0 : planes_with_hit_hor +=1
            if planes_with_hit_hor < 4 : return False
            if planes_with_hit_ver < 4 : return False
        for so in planes:
            for aCl_daughter1 in three_clu[1][so]:
                if aCl_daughter1 in three_clu[2][so] : return False #nOverlap['12']+=1
                aCl_daughter1.GetPosition(A,B)
                if so%2==1:
                    pos_1 = (A[0]+B[0])/2
                else:
                    pos_1 = (A[1]+B[1])/2
                for aCl_daughter2 in three_clu[2][so]:
                    if aCl_daughter2 in three_clu[1][so] :return False
                    if aCl_daughter1==aCl_daughter2 :return False
                    aCl_daughter2.GetPosition(A,B)
                    if so%2==1:
                        pos_2 = (A[0]+B[0])/2
                    else:
                        pos_2 = (A[1]+B[1])/2
#                    if abs(pos_1-pos_2)<0.2: return False

            for aCl_mother in three_clu[0][so]:
                aCl_mother.GetPosition(A,B)
                if aCl_mother in three_clu[1][so] : return False #nOverlap['13']+=1 # return False
                if so%2==1:
                    pos_0 = (A[0]+B[0])/2
                else:
                    pos_0 = (A[1]+B[1])/2
                for aCl_daughter1 in three_clu[1][so]:
                    if aCl_daughter1 in three_clu[0][so] :return False
                    if aCl_daughter1==aCl_mother : return False
                    aCl_daughter1.GetPosition(A,B)
                    if so%2==1:
                        pos_1 = (A[0]+B[0])/2
                    else:
                        pos_1 = (A[1]+B[1])/2
#                    if abs(pos_1-pos_0)<0.2: return False

            for aCl_mother in three_clu[0][so]:
                aCl_mother.GetPosition(A,B)
                if aCl_mother in three_clu[2][so]: return False  #nOverlap['23']+=1 #return False
                if so%2==1:
                    pos_0 = (A[0]+B[0])/2
                else:
                    pos_0 = (A[1]+B[1])/2 
                for aCl_daughter2 in three_clu[2][so]:
                    if aCl_daughter2 in three_clu[0][so] :return False
                    if aCl_daughter2==aCl_mother:return False
                    aCl_daughter2.GetPosition(A,B)
                    if so%2==1:
                        pos_2 = (A[0]+B[0])/2
                    else:
                        pos_2 = (A[1]+B[1])/2
#                    if abs(pos_0-pos_2)<0.2: return False
#        if nOverlap['13']>1 and nOverlap['23']>1 and nOverlap['12']>1 : return False
        return True


 

    def map2Dict(self,hit2point,detID):
        key  = ROOT.std.vector('int')()
        value = ROOT.std.vector('float')()
        theDict = {}
        key.clear()
        value.clear()
        ROOT.fixRoot(hit2point, key, value, detID)
        for k in range(key.size()):
            theDict[key[k]] = value[k]
        return theDict


    def calculate_angle(self,t1,t2):
        v1 = ROOT.TVector3(t1.GetPx(),t1.GetPy(),t1.GetPz())
        v2 = ROOT.TVector3(t2.GetPx(),t2.GetPy(),t2.GetPz())
        cos_theta = v1.Dot(v2)/(t1.GetP()*t2.GetP())
        theta = ROOT.TMath.ACos(cos_theta)
        return theta


    def write_to_txt(self,my_list):
        file_path = "eventlist.txt"
        with open(file_path, "w") as file:
            for item in my_list:
                file.write(f"{item}\n")


    def invariant_mass(self,p1,p2):
        V0Mom=p1+p2
        return V0Mom.M()

    def vertex_array(self,mctrack_array):
        mu_dict={'MuonToMuonPair':[],'GammaToMuonPair':[],'PositronToMuonPair':[]}
        for mctrack in mctrack_array:
            moId = mctrack.GetMotherId()
            if moId<0:continue
            if abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 13 and mctrack.GetProcName() == 'Lepton pair production': mu_dict['MuonToMuonPair'].append(mctrack)
            if abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 22 and mctrack.GetProcName() == 'Lepton pair production': mu_dict['GammaToMuonPair'].append(mctrack)
            if abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 11 and mctrack.GetProcName() == 'Positron annihilation': mu_dict['PositronToMuonPair'].append(mctrack)
        all_empty = all(not mu_dict[key] for key in mu_dict)
        if all_empty : return []
        vertices = {}
        list_of_vertices =[]
        for proc in mu_dict:
            for trident_mu in mu_dict[proc]:
                vx_z = trident_mu.GetStartZ()
                if not vx_z in vertices : vertices[vx_z]=[],proc
                vertices[vx_z][0].append(trident_mu)
        for key in vertices:
            vx_x, vx_y, vx_z = vertices[key][0][0].GetStartX(), vertices[key][0][0].GetStartY(), vertices[key][0][0].GetStartZ()#vertices[key][0][0].GetStartVertex(self.vertex_pos)
            list_of_vertices.append(TridentVertex(vertices[key][1], [vx_x,vx_y,vx_z], vertices[key][0]))
        return list_of_vertices



class TridentVertex:
    def __init__(self, trident_type, vertex_pos, daughters):
        self.trident_type = trident_type
        self.vertex_pos = vertex_pos
        self.daughters = daughters
    def get_trident_type(self):
        return self.trident_type
    def get_vertex_pos(self):
        return self.vertex_pos
    def get_daughters(self):
        return self.daughters


class BaseCut:
    def __init__(self):
        print("snd cuts initiated")

class CutScifiClustersStation(BaseCut):

    def __init__(self):
#        self.threshold_dict = {s * 10 + o: [0, 8] for s in range(1, 6) for o in range(2)}
        self.tools = Tools()

    def pass_cut(self):
        clusters = self.tools.scifi_clusters()
        sorted_clusters = self.tools.scifi_clusters_per_plane(clusters)
        inter_cl_dist = self.tools.inter_cluster_dist(sorted_clusters)
        for key in sorted_clusters:
            ncl = len(sorted_clusters[key])
            if ncl >= 8 : return False
            if key > 31 and ncl < 2 : return False
        if inter_cl_dist[40]<inter_cl_dist[50] or inter_cl_dist[41]<inter_cl_dist[51]: return True
        if inter_cl_dist[40]>inter_cl_dist[50] and inter_cl_dist[41]>inter_cl_dist[51]:return False
        return True

class Tools:
    def __init__(self):
        print("scifi tools initiated")
        
    def scifi_clusters(self):
        trackTask.clusScifi.Clear()
        trackTask.scifiCluster()
        clusters = trackTask.clusScifi
        return clusters

    def scifi_clusters_per_plane(self,clusters):
        sortedClusters = {s * 10 + o: [] for s in range(1, 6) for o in range(2)}
        for aCl in clusters:
            so = aCl.GetFirst()//100000
            sortedClusters[so].append(aCl)
        return sortedClusters

    def inter_cluster_dist(self,clusters_per_plane):
        A,B = ROOT.TVector3(), ROOT.TVector3()
        cluster_seperation_per_plane = {}
        for so in clusters_per_plane:
            cluster_pos = []
            cluster_list = clusters_per_plane[so]
            for aCl in cluster_list:
                aCl.GetPosition(A,B)
                if so%2==1: pos = (A[0]+B[0])/2
                else: pos = (A[1]+B[1])/2
                cluster_pos.append(pos)
            sorted_pos = sorted(cluster_pos,reverse = True)
            if len(sorted_pos)>0 : cluster_seperation_per_plane[so] = abs(sorted_pos[0]-sorted_pos[-1])
        return cluster_seperation_per_plane



x = MuonTridentEngine()
rc = x.efficiency()


