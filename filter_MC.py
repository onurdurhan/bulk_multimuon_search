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

outFile = "filtered_clusters.root"
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




branches_to_copy = {"MCTrack":"ShipMCTrack","ScifiPoint":"ScifiPoints","MuFilterPoint":"MuFilterPoints","EventHeader.":"sndEventHeader","Digi_ScifiHits":"DigiScifiHit_det","Digi_ScifiHits2MCPoints":"DigiScifiHits2MCPoints_det","Cluster_Scifi":"ScifiCluster_det","Digi_MuFilterHits":".DigiMuFilterHit_det","Digi_MuFilterHits2MCPoints":"DigiMuFilterHits2MCPoints_det"}

for b_name in branches_to_copy:
     b=ioman.GetObject(b_name)
     print(b_name,branches_to_copy[b_name])
     ioman.Register(b_name,branches_to_copy[b_name],b,ROOT.kTRUE)
B = ROOT.TList()
B.SetName('BranchList')
for aBranch in branches_to_copy:
    B.Add(ROOT.TObjString(aBranch))
ioman.SetBranchNameList(B)
ioman.WriteFolder()

eventTree.GetEvent(0)
if eventTree.EventHeader.ClassName() == 'SNDLHCEventHeader':
   geo.modules['Scifi'].InitEvent(eventTree.EventHeader)
   geo.modules['MuFilter'].InitEvent(eventTree.EventHeader)
# if faireventheader, rely on user to select correct geofile.


class MuonTridentEngine:
    "generate histograms for Muon Filter"
    def __init__(self):
        print("Initializing Muon Filter Engine ...")
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
        for proc in proc_names:
            weighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0}
            unweighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0}
        p1, p2 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        check_target = lambda aVertex : all(b[0] < p < b[1] for p, b in zip(aVertex, [(-47.7758,-5.5467),(13.0335 ,55.2793),(-1000000,288.)])) ###whole wall : 353.3114
        print("Total nbr of events in this file, ",eventTree.GetEntries())
        for n in range(ioman.GetInTree().GetEntries()):
            if n%1000==0: print('event at %d' % n)
            ioman.GetInTree().GetEvent(n)
            vertices = self.vertex_array(eventTree.MCTrack) 
            w = eventTree.MCTrack[0].GetWeight()
            h2p = eventTree.Digi_ScifiHits2MCPoints.First()
            multiple_entries = 0
            clusters = self.tools.scifi_clusters()
            recoble_event = False
            for aVx in vertices:
                if aVx.get_vertex_pos()[2]>280. : continue 
                daughters = aVx.get_daughters()
                multiple_entries+=1
                unweighted[aVx.get_trident_type()]['norm']+=1
                weighted[aVx.get_trident_type()]['norm']+=w
                if len(daughters)!=2 : continue
                three_clu = self.three_clusters(h2p, daughters[0], daughters[1], clusters)
                if self.check_three_accepted_clusterwise(three_clu):
                    unweighted[aVx.get_trident_type()]['inAcceptance']+=1
                    weighted[aVx.get_trident_type()]['inAcceptance']+=w
                    if self.with_original_muon(three_clu):
                        unweighted[aVx.get_trident_type()]['with_original_muon']+=1
                        weighted[aVx.get_trident_type()]['with_original_muon']+=w
                    else : continue
                    if self.check_all_planes(clusters): #self.check_planes_with_hits() :#self.check_all_planes(clusters):
                        unweighted[aVx.get_trident_type()]['hits_exist_in_all']+=1
                        weighted[aVx.get_trident_type()]['hits_exist_in_all']+=w 
                        if self.recoble(three_clu): recoble_event = True
            if recoble_event : ioman.Fill()
        ioman.Write()
        return True

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


