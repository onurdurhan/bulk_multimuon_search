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
import cppyy
from scipy.optimize import fsolve
from scipy.optimize import minimize
from itertools import combinations
from Tools import Tools
import BaseCut as BaseCut

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file data and MC",default="",required=False)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geofile", default=os.environ["EOSSHIP"]+"/eos/experiment/sndlhc/users/odurhan/MuonicTridentMC/geofile_full.Ntuple-TGeant4_boost100.0.root")
parser.add_argument("-j","--jobID",dest="jobID",help="job ID",default=0,required=False)
parser.add_argument("-t","--trackType",dest="trackType",help="scifi of ds",default="scifi",required=True)
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


geo_file =options.geoFile

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
        ut.bookHist(self.hist,"xz","",200,-7500.,2500.,200,-100.,100.)
        ut.bookHist(self.hist,"yz","",200,-7500.,2500.,200, -100.,100.)
        ut.bookHist(self.hist,"xy","",200,-100,100,200,-100,100)
        for o in range(2):
            for s in range(1,6):
                ut.bookHist(self.hist,"n_clu_so"+str(s*10+o)+"_MC","scifi hit density "+ str(s*10+o),40,0,40)
            for s in range(1,5):
                ut.bookHist(self.hist,'n_hitsDS_so'+str(s*10+o)+"_MC",'ds hit density' + str(s*10+o),10,0,10)
        self.p1, self.p2 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        for key in self.trident_type:
            ut.bookHist(self.hist,"inv_"+key,"",25,0.,1.)
        self.cases = {1:"single",2:"dimuon",3:"threemuon"}
        self.materials = {1:'Rock',2:'Concrete',3:'air'}
        ut.bookHist(self.hist,'charge','',29,-14,14)
        dmin_bins = [0,1.,2.5,4.5,7.5,11.5,16.5,25.]  # Variable bin widths in Y
        # Number of bins = (number of edges - 1)
        n_dmin_bins = len(dmin_bins) - 1
        dmin_bins_array = ROOT.std.vector('double')(dmin_bins)
        self.dmax_vs_dmin={}
        self.z_ranges=[7150,5000,3000,1000,0]
        self.E_ranges=[30,500,1000,3000,5000]
        for key in self.trident_type:
            ut.bookHist(self.hist,'z_total'+key,'',20,-7500.,280.)
            ut.bookHist(self.hist,'z_pass'+key,'',20,-7500.,280)
            ut.bookHist(self.hist,'z_total_acc'+key,'',20,-7500,500)
            ut.bookHist(self.hist,'z_pass_acc'+key,'',20,-7500,500)
            ut.bookHist(self.hist,'E_pass_acc'+key,'',20,0,5000)
            ut.bookHist(self.hist,'E_total_acc'+key,'',20,0,5000)
            ut.bookHist(self.hist,'E_total'+key,'',20,0,5000)
            ut.bookHist(self.hist,'E_pass'+key,'',20,0,5000)
            ut.bookHist(self.hist,'material_'+key,'',3,0,3)
            ut.bookHist(self.hist,"E_daughters_"+key,'',100,0,3000)
            ut.bookHist(self.hist,f"theta_for_ship_{key}","",100,0,100)
            if key =='MuonToMuonPair': 
                max_mass = 1.
            else :
                max_mass = 10.
            ut.bookHist(self.hist,f"inv_mass_for_ship_{key}","",100,0,max_mass)

            for b in range(1,4):
                self.hist['material_'+key].GetXaxis().SetBinLabel(b,self.materials[b])
            if key!='MuonToMuonPair':
                ut.bookHist(self.hist,'theta_total'+key,'',10,0,50)
                ut.bookHist(self.hist,'theta_pass' +key,'',10,0,50)
                ut.bookHist(self.hist,'E_vs_theta_pass'+key,'',25,0,50,25,0,5000)
                ut.bookHist(self.hist,'E_vs_theta_total'+key,'',25,0,50,25,0,5000)
                ut.bookHist(self.hist,'z_vs_theta_pass'+key,'',25,0,25,25,-8000,0)
                ut.bookHist(self.hist,'z_vs_theta_total'+key,'',25,0,25,25,-8000,0)
            else:
                ut.bookHist(self.hist,'theta_total'+key,'',10,0,5)
                ut.bookHist(self.hist,'theta_pass'+key,'',10,0,5)
                ut.bookHist(self.hist,'E_vs_theta_pass'+key,'',50,0,5,50,0,5000)
                ut.bookHist(self.hist,'E_vs_theta_total'+key,'',50,0,5,50,0,5000)
                ut.bookHist(self.hist,'z_vs_theta_pass'+key, '',50,0,5,50,-8000,0)
                ut.bookHist(self.hist,'z_vs_theta_total'+key,'',50,0,5,50,-8000,0)
            for p in range(2):
                if key == 'MuonToMuonPair': 
                    nbin=25
                else : nbin = 25
                for i in range(len(self.z_ranges)-1):
                    self.dmax_vs_dmin[key+'_histo_dmax_vs_dmin_'+str(p)+"_"+str(self.z_ranges[i])+"_"+str(self.z_ranges[i+1])] = ROOT.TH2D(key+"h2"+str(p)+"_"+str(self.z_ranges[i])+"_"+str(self.z_ranges[i+1]), "2D Histogram with Variable Y Binning;X axis;Y axis",n_dmin_bins,dmin_bins_array.data(),nbin,0,25.)
                    self.dmax_vs_dmin[key+'_histo_dmax_vs_dmin_'+str(p)+"_"+str(self.E_ranges[i])+"_"+str(self.E_ranges[i+1])] = ROOT.TH2D(key+"h2"+str(p)+"_"+str(self.E_ranges[i])+"_"+str(self.E_ranges[i+1]), "2D Histogram with Variable Y Binning;X axis;Y axis",n_dmin_bins,dmin_bins_array.data(),nbin,0,25.)
                    ut.bookHist(self.hist,f'{key}_dmin_dmax_ratio_{p}_{self.z_ranges[i]}_{self.z_ranges[i+1]}',"dmin to dmax ratio",60,0,0.6)
                    ut.bookHist(self.hist,f'{key}_dmin_dmax_ratio_{p}_{self.E_ranges[i]}_{self.E_ranges[i+1]}',"dmin to dmax ratio",30,0,0.6)
 
        ut.bookHist(self.hist,"E_loss","",50,0,20)

        
        self.scifi_hit_density_cut= BaseCut.ScifiHitDensityCut(eventTree,geo)
        self.scifi_fiducial_cut = BaseCut.FiducialCut()
        self.tools = Tools(trackTask,geo)        

    def efficiency(self,trackType):
        proc_names = ['MuonToMuonPair','GammaToMuonPair','PositronToMuonPair']
        weighted = {}
        unweighted = {}
        cases = {1:"single",2:"dimuon",3:"threemuon"}
        for proc in proc_names:
            weighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0, 'fiducial':0}
            unweighted[proc]={'norm':0,'inAcceptance':0,'with_original_muon':0,'hits_exist_in_all':0,'recoble':0,'3/3':0,'2/3':0,'2/3':0,'2/2':0,'sparse':0,'fiducial':0}
        print("Total nbr of events in this file, ",eventTree.GetEntries())
        A,B = ROOT.TVector3(),ROOT.TVector3()
        multiple_entries = 0
        single_muon_tracks = 0
        not_convergent = 0
        events_with_reco_tracks=[]
        for n in range(eventTree.GetEntries()):
            if n%10000==0: print('event at %d' % n)
            eventTree.GetEvent(n)
            vertices = self.vertex_array(eventTree.MCTrack)
            w = eventTree.MCTrack[0].GetWeight()
            """
            for single muon tracks
            """
            for aTrack in Reco_MuonTracks: aTrack.Delete()
            Reco_MuonTracks.Clear("C")
            rc = trackTask.ExecuteTask("Scifi")
            if not Reco_MuonTracks.GetEntries()==1: 
                continue
            theTrack = Reco_MuonTracks[0]
            if not hasattr(theTrack,"getFittedState"):continue
            status = theTrack.getFitStatus()
            if not status.isFitConverged() :
                not_convergent+=w
                continue
            zPos = self.tools.getAverageZpositions()
            testPlane=51
            z = zPos['Scifi'][testPlane]
            parallelToZ = ROOT.TVector3(0,0,1.)
            rep     = ROOT.genfit.RKTrackRep(13)
            state  = ROOT.genfit.StateOnPlane(rep)
            # find closest track state
            mClose = 0
            mZmin = 999999.
            for m in range(0,theTrack.getNumPointsWithMeasurement()):
                st   = theTrack.getFittedState(m)
                Pos = st.getPos()
                if abs(z-Pos.z())<mZmin:
                    mZmin = abs(z-Pos.z())
                    mClose = m
            if mZmin>10000:
                print("something wrong here with measurements",mClose,mZmin,theTrack.getNumPointsWithMeasurement())
                continue
            fstate =  theTrack.getFittedState(mClose)
            pos,mom = fstate.getPos(),fstate.getMom()
            rep.setPosMom(state,pos,mom)
            NewPosition = ROOT.TVector3(0., 0., z)   # assumes that plane in global coordinates is perpendicular to z-axis, which is not true for TI18 geometry.
            rep.extrapolateToPlane(state, NewPosition, parallelToZ )
            pos = state.getPos()
            xEx,yEx = pos.x(),pos.y()
            Simona_vol = True
            if yEx<=18. or yEx >= 49:
                Simona_vol = False
            if xEx<=-42 or xEx>=-11:
                Simona_vol = False
            if Simona_vol: 
                single_muon_tracks+=w
            events_with_reco_tracks.append(n)
            continue
            if trackType == "scifi":
                h2p = eventTree.Digi_ScifiHits2MCPoints.First()
            else :
                h2p = eventTree.Digi_MuFilterHits2MCPoints.First()
            if trackType == "scifi":
                clusters =self.tools.scifi_clusters()
            recoble_event = False
            tracks_counted = False
            n3D = [0,0]
            if len(vertices)>1:multiple_entries+=w

            for aVx in vertices:
                x,y,z =  aVx.get_vertex_pos()[0], aVx.get_vertex_pos()[1],aVx.get_vertex_pos()[2]
                """
                for ship
                if eventTree.MCTrack[0].GetEnergy()>0.:
                    p1,p2 = ROOT.TLorentzVector(),ROOT.TLorentzVector()
                    theta_ship = self.calculate_angle(aVx.get_daughters()[0],aVx.get_daughters()[1])
                    self.hist[f"theta_for_ship_{aVx.get_trident_type()}"].Fill(theta_ship*1E3)
                    found = False
                    for muon in aVx.get_daughters():
                        muon.Get4Momentum(p1)
                        m =105.66/1E3
                        for p in eventTree.ScifiPoint:
                            if p.GetTrackID()<0:continue
                            px = p.GetPx()
                            py = p.GetPy()
                            pz = p.GetPz()
                            E_f = (px**2 + py**2 + pz**2 + m**2)**0.5
                            if eventTree.MCTrack[p.GetTrackID()] == aVx.get_mother():
                                p3 = ROOT.TLorentzVector(px, py, pz, E_f)
                                found = True
                            if found : break

                    mass=self.invariant_mass(p1,p2,p3)
                    self.hist[f"inv_mass_for_ship_{aVx.get_trident_type()}"].Fill(mass)
                    print(mass, aVx.get_trident_type(),aVx.get_mother().GetEnergy())
                """    

                if z>280:continue
                if aVx.get_trident_type() == 'MuonToMuonPair':
                    self.hist["xy"].Fill(x,y,w)
                    self.hist["yz"].Fill(z,y,w)
                    self.hist["xz"].Fill(z,x,w)
                daughters = aVx.get_daughters()
                for daughter in daughters:
                    self.hist['E_daughters_'+aVx.get_trident_type()].Fill(daughter.GetEnergy(),w)
                unweighted[aVx.get_trident_type()]['norm']+=1
                weighted[aVx.get_trident_type()]['norm']+=w
                if len(daughters)!=2 : 
                    continue
                for p in eventTree.ScifiPoint:
                    if p.GetTrackID()<0:continue
                    mc_track=eventTree.MCTrack[p.GetTrackID()]
#                    if not mc_track.GetMotherId()== 0: continue
                    if mc_track in daughters and mc_track.GetStartZ()>-1100. and mc_track.GetStartZ()<-999.:
                        px = p.GetPx()
                        py = p.GetPy()
                        pz = p.GetPz()    
                        m=105.66/1000
                        E_f = (px**2 + py**2 + pz**2 + m**2)**0.5
                        self.hist["E_loss"].Fill(mc_track.GetEnergy()-E_f)
                if trackType == "scifi":
                    three_clu = self.three_clusters(h2p, daughters[0], daughters[1], clusters,trackType)
                else :
                    three_clu = self.three_clusters(h2p,daughters[0],daughters[1],eventTree.Digi_MuFilterHits,trackType)
                inAcc = self.check_three_accepted(three_clu)
                self.hist['z_total_acc'+aVx.get_trident_type()].Fill(aVx.get_daughters()[0].GetStartZ(),w)
                self.hist['E_total_acc'+aVx.get_trident_type()].Fill(aVx.get_mother().GetEnergy(),w)
                if not inAcc:
                    continue
                unweighted[aVx.get_trident_type()]['inAcceptance']+=1
                weighted[aVx.get_trident_type()]['inAcceptance']+=w
                if self.with_original_muon(three_clu):
                    unweighted[aVx.get_trident_type()]['with_original_muon']+=1
                    weighted[aVx.get_trident_type()]['with_original_muon']+=w
                if trackType == "scifi": 
                    all_planes = self.check_all_planes(clusters,trackType)
                else : 
                    all_planes = self.check_all_planes(eventTree.Digi_MuFilterHits,trackType)
                if not all_planes: 
                    continue
                unweighted[aVx.get_trident_type()]['hits_exist_in_all']+=1
                weighted[aVx.get_trident_type()]['hits_exist_in_all']+=w 
                if trackType== "scifi":
                    if self.recoble(three_clu): recoble_event = True
                else :
                    if self.recoble_DS(three_clu): recoble_event = True        
                if not recoble_event: continue
                self.hist['z_pass_acc'+aVx.get_trident_type()].Fill(aVx.get_daughters()[0].GetStartZ(),w)
                self.hist['E_pass_acc'+aVx.get_trident_type()].Fill(aVx.get_mother().GetEnergy(),w)
                unweighted[aVx.get_trident_type()]['recoble']+=1
                weighted[aVx.get_trident_type()]['recoble']+=w 
                if tracks_counted : continue   ## count tracks once for all vertices 
                opening_angle = self.calculate_angle(aVx.get_daughters()[0],aVx.get_daughters()[1])
                self.hist['E_total'+aVx.get_trident_type()].Fill(aVx.get_mother().GetEnergy(),w)
                self.hist['z_total'+aVx.get_trident_type()].Fill(aVx.get_daughters()[0].GetStartZ(),w)
                opening_angle = self.calculate_angle(aVx.get_daughters()[0],aVx.get_daughters()[1])
                self.hist['theta_total'+aVx.get_trident_type()].Fill(opening_angle*1000,w)
                self.hist['E_vs_theta_total'+aVx.get_trident_type()].Fill(opening_angle*1000,eventTree.MCTrack[aVx.get_daughters()[0].GetMotherId()].GetEnergy(),w)   
                self.hist['z_vs_theta_total'+aVx.get_trident_type()].Fill(opening_angle*1000,aVx.get_daughters()[0].GetStartZ(),w)
                if trackType == "scifi":
                    trackTask.multipleTrackCandidates(nMaxCl=8,dGap=0.2,dMax=0.8,dMax3=0.8,ovMax=1,doublet=True,debug=False)
                else:
                    trackTask.multipleTrackCandidates2(geo.modules['MuFilter'],eventTree.Digi_MuFilterHits, debug = False)
                reco_2dTracks = {"scifi":{0:[],1:[]}, "ds" :{0:[],1:[]}}
                for p in range(2):
                    for trackId in trackTask.multipleTrackStore['trackCand'][p]:
                        if trackId in trackTask.multipleTrackStore['cloneCand'][p]: continue
                        n3D[p]+=1
                        if trackType == "scifi":
                            reco_2dTracks["scifi"][p].append(trackTask.multipleTrackStore['trackCand'][p][trackId])
                        else :
                            reco_2dTracks["ds"][p].append(trackTask.multipleTrackStore['trackCand'][p][trackId])
                if n3D == [3,3] :
                    unweighted[aVx.get_trident_type()]['3/3']+=1
                    weighted[aVx.get_trident_type()]['3/3']+=w
                if  n3D == [2,3] or n3D == [3,2]:
                    unweighted[aVx.get_trident_type()]['2/3']+=1
                    weighted[aVx.get_trident_type()]['2/3']+=w
                if n3D == [2,2]:
                    unweighted[aVx.get_trident_type()]['2/2']+=1
                    weighted[aVx.get_trident_type()]['2/2']+=w
                tracks_counted = True
                if not n3D == [3,3]: continue ### now focus on 3-track events !        
                if trackType == "scifi":
                    sparse_event = self.scifi_hit_density_cut.pass_cut(reco_2dTracks)
                    if not sparse_event:
                        continue
                    unweighted[aVx.get_trident_type()]['sparse']+=1
                    weighted[aVx.get_trident_type()]['sparse']+=w
                    scifi_fiducial = self.scifi_fiducial_cut.pass_cut(reco_2dTracks)
                    if not scifi_fiducial :
                        continue
                    unweighted[aVx.get_trident_type()]['fiducial']+=1
                    weighted[aVx.get_trident_type()]['fiducial']+=w 
                    self.hist['E_pass'+aVx.get_trident_type()].Fill(eventTree.MCTrack[aVx.get_daughters()[0].GetMotherId()].GetEnergy(),w)
                    self.hist['z_pass'+aVx.get_trident_type()].Fill(aVx.get_daughters()[0].GetStartZ(),w)
                    opening_angle = self.calculate_angle(aVx.get_daughters()[0],aVx.get_daughters()[1])
                    self.hist['theta_pass'+aVx.get_trident_type()].Fill(opening_angle*1000,w)
                    self.hist['E_vs_theta_pass'+aVx.get_trident_type()].Fill(opening_angle*1000,eventTree.MCTrack[aVx.get_daughters()[0].GetMotherId()].GetEnergy(),w) 
                    self.hist['z_vs_theta_pass'+aVx.get_trident_type()].Fill(opening_angle*1000,aVx.get_daughters()[0].GetStartZ(),w)
                    self.hist['charge'].Fill(eventTree.MCTrack[0].GetPdgCode(),w)            
                    thetas=[]
                    distances = []
                    dmins=[]
                    dmax=[]
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
                                X1,X2 = line1.Eval(300),line2.Eval(300)
                                slope1,slope2=line1.GetParameter(1),line2.GetParameter(1)
                                theta=self.angle_between_lines(slope1, slope2)
                                d=abs(X1-X2)
                                if d not in pairs :
                                    pairs[d]=theta
                        sorted_dict = dict(sorted(pairs.items()))
                        first_key = next(iter(sorted_dict))
                        last_key = next(reversed(sorted_dict))
                        dmin_dmax_ratio = first_key/last_key
                        dmins.append(first_key)
                        dmax.append(last_key)
                        for i in range(len(self.z_ranges)-1):
                            if -self.z_ranges[i]<aVx.get_vertex_pos()[2] and -self.z_ranges[i+1]>aVx.get_vertex_pos()[2]:
                                self.dmax_vs_dmin[f'{aVx.get_trident_type()}_histo_dmax_vs_dmin_{p}_{self.z_ranges[i]}_{self.z_ranges[i+1]}'].Fill(first_key,last_key,w)
#                                self.hist[f'{aVx.get_trident_type()}_dmin_dmax_ratio_{p}_{self.z_ranges[i]}_{self.z_ranges[i+1]}'].Fill(dmin_dmax_ratio,w)
                        for i in range(len(self.E_ranges)-1):
                            if self.E_ranges[i]<aVx.get_mother().GetEnergy() and self.E_ranges[i+1]>aVx.get_mother().GetEnergy():
                                self.dmax_vs_dmin[f'{aVx.get_trident_type()}_histo_dmax_vs_dmin_{p}_{self.E_ranges[i]}_{self.E_ranges[i+1]}'].Fill(first_key,last_key,w)
#                                self.hist[f'{aVx.get_trident_type()}_dmin_dmax_ratio_{p}_{self.E_ranges[i]}_{self.E_ranges[i+1]}'].Fill(dmin_dmax_ratio,w)
 
                        vx_x, vx_y, vx_z = daughters[0].GetStartX(),daughters[0].GetStartY(),daughters[0].GetStartZ()
                        node = geo.sGeo.FindNode(vx_x,vx_y,vx_z)
                        material_name = node.GetMedium().GetMaterial().GetName()
                        if material_name == 'Rock' or material_name == 'MolasseRock':
                            self.hist['material_'+aVx.get_trident_type()].Fill(0.5,w)
                        elif material_name == 'Concrete':
                            self.hist['material_'+aVx.get_trident_type()].Fill(1.5,w)
                        elif material_name == 'air':
                            self.hist['material_'+aVx.get_trident_type()].Fill(2.5,w)
                        else:
                            print('different material', material_name)
                    dmin_dmax_ratio=ROOT.TMath.Sqrt(dmins[0]**2+dmins[1]**2)/ROOT.TMath.Sqrt(dmax[0]**2+dmax[1]**2)
                    for i in range(len(self.z_ranges)-1):
                        if -self.z_ranges[i]<aVx.get_vertex_pos()[2] and -self.z_ranges[i+1]>aVx.get_vertex_pos()[2]:
                            self.hist[f'{aVx.get_trident_type()}_dmin_dmax_ratio_{p}_{self.z_ranges[i]}_{self.z_ranges[i+1]}'].Fill(dmin_dmax_ratio,w)

            continue
            if len(vertices)>1 : continue
            for aHit in eventTree.Digi_MuFilterHits:
                system = aHit.GetSystem()
                if system==3 :continue
                detID = aHit.GetDetectorID()
                the_dict=self.map2Dict(eventTree.Digi_MuFilterHits2MCPoints.First(),detID)
                n_3 = []
                for key in the_dict:
                    trackID = eventTree.MuFilterPoint[key].GetTrackID()
                    if abs(eventTree.MCTrack[trackID].GetPdgCode())==13:
                        n_3.append(trackID)
                if len(n_3)<1 or len(n_3)>3:continue
                sumSignal = self.map2Dict_qdc(aHit,'SumOfSignals')
                plane  = (detID%10000)//1000
                bar = detID%1000
                if system == 1 :
                    bar_key = "Vetoplane"+str(plane)+"horizontalbar"+str(bar)
                if system == 2:
                    bar_key = "USplane"+str(plane)+"horizontalbar"+str(bar)
#                self.hist["qdc_"+cases[len(n_3)]+bar_key].Fill(sumSignal['Sum'])
                
        print("UNWEIGHTED ", unweighted)
        print("WEIGHTED: ", weighted)
        print("Multple entries ", multiple_entries)
        print("single muon tracks ", single_muon_tracks)
        print("non-convergent tracks",not_convergent)
        self.write_to_txt(events_with_reco_tracks)
        # Write dictionary to a file in JSON format
        with open('data_unweighted_'+str(options.jobID)+'.json', 'w') as json_file:
            json.dump(unweighted, json_file)
        with open('data_weighted_'+str(options.jobID)+'.json','w') as json_file:
            json.dump(weighted, json_file)

#        with open('data'+str(options.jobID)+'.json', 'r') as json_file:
#            loaded_dict = json.load(json_file)
#            print(loaded_dict)
        ut.writeHists(self.hist,"histos_newdigi"+str(options.jobID)+"_"+options.trackType+".root")
        ana_file=ROOT.TFile("histos_newdigi"+str(options.jobID)+"_"+options.trackType+".root","UPDATE")
        ana_file.cd()
        """
        for p in range(2):
            for proc in proc_names:
                graph_name=proc+'_dmax_vs_dmin_'+str(p)
                self.gr[graph_name]=ROOT.TGraph(len(dmax[proc+"_"+str(p)]),dmin[proc+"_"+str(p)],dmax[proc+"_"+str(p)])
        """
        for p in range(2):
            for proc in proc_names:
#                graph_name=proc+'_dmax_vs_dmin_'+str(p)
#                self.gr[graph_name].Write(graph_name)
                for i in range(len(self.z_ranges)-1):
                    self.dmax_vs_dmin[proc+"_histo_dmax_vs_dmin_"+str(p)+"_"+str(self.z_ranges[i])+"_"+str(self.z_ranges[i+1])].Write(proc+"_histo_dmax_vs_dmin_"+str(p)+"_"+str(self.z_ranges[i])+"_"+str(self.z_ranges[i+1]))
                for i in range(len(self.E_ranges)-1):
                    self.dmax_vs_dmin[proc+"_histo_dmax_vs_dmin_"+str(p)+"_"+str(self.E_ranges[i])+"_"+str(self.E_ranges[i+1])].Write(proc+"_histo_dmax_vs_dmin_"+str(p)+"_"+str(self.E_ranges[i])+"_"+str(self.E_ranges[i+1]))
 
        ana_file.Close()








    def Q2(self,vertex):
        """
        mother = eventTree.MCTrack[vertex.get_daughters()[0].GetMotherId()]
        P_i = ROOT.TLorentzVector()
        mother.Get4Momentum(P_i)
        m = 105.66/1000
        for p in eventTree.ScifiPoint:
            if p.GetTrackID()==0:
                px = p.GetPx()
                py = p.GetPy()
                pz = p.GetPz()
                E_f = (px**2 + py**2 + pz**2 + m**2)**0.5
                break
        P_f = ROOT.TLorentzVector(px, py, pz, E_f)
        Q = P_i-P_f
        Q2 = -Q.M2()
        return Q.Energy()
        return Q2
        """
        P1, P2 = ROOT.TLorentzVector(), ROOT.TLorentzVector()
        vertex.get_daughters()[0].Get4Momentum(P1)
        vertex.get_daughters()[1].Get4Momentum(P2)
        Q = P1+P2
        Q2 = Q.M2()
        Q2 = Q.Energy()
        return Q2

    def calculate_opening_angles_projections(self,p_mu_plus, p_mu_minus):
        """
        Calculate the opening angles between muon pair in the xz and yz planes.
    
        Parameters:
        p_mu_plus  -- Momentum vector of mu+ as a tuple (px, py, pz)
        p_mu_minus -- Momentum vector of mu- as a tuple (px, py, pz)
    
        Returns:
        theta_xz -- Opening angle in the xz plane (in radians)
        theta_yz -- Opening angle in the yz plane (in radians)
        """
        # Extract components of the momentum vectors
        p_plus_x, p_plus_y, p_plus_z = p_mu_plus
        p_minus_x, p_minus_y, p_minus_z = p_mu_minus

        # Calculate the dot products and magnitudes for xz plane
        dot_product_xz = p_plus_x * p_minus_x + p_plus_z * p_minus_z
        magnitude_plus_xz = np.sqrt(p_plus_x**2 + p_plus_z**2)
        magnitude_minus_xz = np.sqrt(p_minus_x**2 + p_minus_z**2)

        # Calculate the opening angle in the xz plane
        theta_xz = np.arccos(dot_product_xz / (magnitude_plus_xz * magnitude_minus_xz))

        # Calculate the dot products and magnitudes for yz plane
        dot_product_yz = p_plus_y * p_minus_y + p_plus_z * p_minus_z
        magnitude_plus_yz = np.sqrt(p_plus_y**2 + p_plus_z**2)
        magnitude_minus_yz = np.sqrt(p_minus_y**2 + p_minus_z**2)

        # Calculate the opening angle in the yz plane
        theta_yz = np.arccos(dot_product_yz / (magnitude_plus_yz * magnitude_minus_yz))

        return theta_xz, theta_yz



    def angle_between_lines(self,slope1, slope2):
        if slope1 == slope2:
            return 0.0  # The lines are parallel
        # Calculate the tangent of the angle between the lines
        tan_theta = abs((slope2 - slope1) / (1 + slope1 * slope2))
    
        # Convert the tangent to an angle in degrees
        angle = ROOT.TMath.ATan(tan_theta) #* (180 / ROOT.TMath.Pi())
    
        return angle


    def is_UShit_isolated(self,theHit):
        detID = theHit.GetDetectorID()
        the_plane  = (detID%10000)//1000
        the_bar = detID%1000
        isolated = True
        for aHit in eventTree.Digi_MuFilterHits:
            if aHit.GetSystem()!=2 : continue
            if aHit.GetDetectorID() == detID : continue 
            if (aHit.GetDetectorID()%10000)//1000 == the_plane:
                if abs(the_bar-aHit.GetDetectorID()%1000)==1:isolated = False
        return isolated


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



    def three_clusters(self,h2p,daughter1,daughter2,obj_array,trackType):
        if trackType=='scifi':
            sorted_clusters = self.tools.scifi_clusters_per_plane(obj_array)
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
        else:
            sorted_hits = {s * 10 + o: [] for s in range(1, 5) for o in range(2)}
            for aHit in obj_array:
                if aHit.GetSystem() !=3 : continue
                detID = aHit.GetDetectorID()
                station  = (detID%10000)//1000  # station number
                bar = (detID%1000)
                so = (station+1)*10+(bar>59)
                sorted_hits[so].append(aHit)
            three_clu = {k: {s * 10 + o: [] for s in range(1, 5) for o in range(2)} for k in range(3)} 
            for so in sorted_hits:
                for aHit in sorted_hits[so]:
                    detID = aHit.GetDetectorID()
                    theDict = self.map2Dict(h2p,detID)
                    for link_point in theDict:
                        if eventTree.MuFilterPoint[link_point].GetTrackID()< -1 :continue
                        registered_track = eventTree.MCTrack[eventTree.MuFilterPoint[link_point].GetTrackID()]
                        if registered_track == daughter1 and aHit not in three_clu[1][so] : three_clu[1][so].append(aHit)
                        if registered_track == daughter2 and aHit not in three_clu[2][so] : three_clu[2][so].append(aHit)
                        if registered_track == eventTree.MCTrack[0] and aHit not in three_clu[0][so] : three_clu[0][so].append(aHit)
            return three_clu 

    def check_three_accepted(self, three_clu):
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




    def check_all_planes(self,obj_array,trackType):
        if trackType=='scifi':
            sorted_clusters = self.tools.scifi_clusters_per_plane(obj_array)
            for so in sorted_clusters:
                if len(sorted_clusters[so]) < 1 :return False
            return True

        else:
            planes_with_hits = {s*10+o:False for s in range(1,5) for o in range(2)}
            for aHit in obj_array:
                if aHit.GetSystem()!=3 : continue
                detID =aHit.GetDetectorID()
                station  = (detID%10000)//1000  # station number
                bar = (detID%1000)
                so = (station+1)*10+(bar>59)
                planes_with_hits[so]=True
            planes_with_hits[40]=True
            return all(planes_with_hits.values())



    def recoble_DS(self,three_clu):
        planes = [s * 10 + o for s in range(1, 5) for o in range(2)]
        nOverlap = {'12':0,'13':0,'23':0}
        for k in three_clu:
            planes_with_hit_hor = 0
            planes_with_hit_ver = 0
            for so in three_clu[k]:
                if so%2==1:
                    if len(three_clu[k][so])>0 : planes_with_hit_ver += 1
                else:
                    if len(three_clu[k][so])>0 : planes_with_hit_hor +=1
            if planes_with_hit_hor < 3 : return False
            if planes_with_hit_ver < 4 : return False
        for so in planes:
            if so == 40 : continue
            for aHit_daughter1 in three_clu[1][so]:
                if aHit_daughter1 in three_clu[2][so] : return False #nOverlap['12']+=1
                for aHit_daughter2 in three_clu[2][so]:
                    if aHit_daughter2 in three_clu[1][so] :return False
                    if aHit_daughter1==aHit_daughter2 :return False

            for aHit_mother in three_clu[0][so]:
                if aHit_mother in three_clu[1][so] : return False #nOverlap['13']+=1 # return False
                for aHit_daughter1 in three_clu[1][so]:
                    if aHit_daughter1 in three_clu[0][so] :return False
                    if aHit_daughter1==aHit_mother : return False

            for aHit_mother in three_clu[0][so]:
                if aHit_mother in three_clu[2][so]: return False  #nOverlap['23']+=1 #return False
                for aHit_daughter2 in three_clu[2][so]:
                    if aHit_daughter2 in three_clu[0][so] :return False
                    if aHit_daughter2==aHit_mother:return False

#        if nOverlap['13']>1 and nOverlap['23']>1 and nOverlap['12']>1 : return False
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
                for aCl_daughter2 in three_clu[2][so]:
                    if aCl_daughter2 in three_clu[1][so] :return False
                    if aCl_daughter1==aCl_daughter2 :return False
            for aCl_mother in three_clu[0][so]:
                aCl_mother.GetPosition(A,B)
                if aCl_mother in three_clu[1][so] : return False #nOverlap['13']+=1 # return False
                for aCl_daughter1 in three_clu[1][so]:
                    if aCl_daughter1 in three_clu[0][so] :return False
                    if aCl_daughter1==aCl_mother : return False
            for aCl_mother in three_clu[0][so]:
                aCl_mother.GetPosition(A,B)
                if aCl_mother in three_clu[2][so]: return False  #nOverlap['23']+=1 #return False
                for aCl_daughter2 in three_clu[2][so]:
                    if aCl_daughter2 in three_clu[0][so] :return False
                    if aCl_daughter2==aCl_mother:return False
#                    if abs(pos_0-pos_2)<0.2: return False
#        if nOverlap['13']>1 and nOverlap['23']>1 and nOverlap['12']>1 : return False
        return True

    def with_original_muon(self,three_clu):
        n_original_mu = 0
        for so in three_clu[0]:
            n_original_mu+=len(three_clu[0][so])
        if n_original_mu<1: return False
        else : return True




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


    def invariant_mass(self,p1,p2,p3):
        V0Mom=p1+p2+p3
        return V0Mom.M()

    def vertex_array(self,mctrack_array):
        mu_dict={'MuonToMuonPair':[],'GammaToMuonPair':[],'PositronToMuonPair':[]}
        for mctrack in mctrack_array:
            moId = mctrack.GetMotherId()
            if moId<0:
                continue
            if   abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 13 and mctrack.GetProcName() == 'Lepton pair production':
                mu_dict['MuonToMuonPair'].append(mctrack)
            elif abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 22 and mctrack.GetProcName() == 'Lepton pair production': 
                mu_dict['GammaToMuonPair'].append(mctrack)
            elif abs(mctrack.GetPdgCode())==13 and abs(mctrack_array[moId].GetPdgCode()) == 11 and mctrack.GetProcName() == 'Positron annihilation':
                mu_dict['PositronToMuonPair'].append(mctrack)
            else:
                continue
        all_empty = all(not mu_dict[key] for key in mu_dict)
        if all_empty :
            return []
        list_of_vertices =[]
        for proc in mu_dict:
            keys = list(mu_dict[proc])
            for i in range(len(mu_dict[proc])-1):
                for j in range(i+1,len(mu_dict[proc])):
                    mu1,mu2=keys[i],keys[j]
                    moId_1, moId_2= mu1.GetMotherId(), mu2.GetMotherId()
                    z1,z2 = mu1.GetStartZ(), mu2.GetStartZ()
                    if moId_1 == moId_2 and abs(z1-z2)<1e-6:
                        vx_x, vx_y, vx_z = mu1.GetStartX(), mu1.GetStartY(), mu1.GetStartZ()
                        list_of_vertices.append(TridentVertex(proc, [vx_x,vx_y,vx_z], [mu1,mu2], mctrack_array[moId_1]))
        return list_of_vertices



class TridentVertex:
    def __init__(self, trident_type, vertex_pos, daughters,mother):
        self.trident_type = trident_type
        self.vertex_pos = vertex_pos
        self.daughters = daughters
        self.mother = mother
    def get_trident_type(self):
        return self.trident_type
    def get_vertex_pos(self):
        return self.vertex_pos
    def get_daughters(self):
        return self.daughters
    def get_mother(self):
        return self.mother



x = MuonTridentEngine()
rc = x.efficiency(options.trackType)
