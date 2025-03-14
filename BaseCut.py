import ROOT
class BaseCut:
    def __init__(self):        
        print("snd cuts initiated")

class ScifiHitDensityCut(BaseCut):
    def __init__(self,tree,geo,rho=10):
        print("Cut for number of  hits near the track")
        self.rho_cut=rho
        self.eventTree = tree
        self.geo = geo
 
    def pass_cut(self, reco_2dTracks):
        passes = True
        A,B = ROOT.TVector3(), ROOT.TVector3()
        noisy_stations = {s:0 for s in range(1,6)}
        for p in range(2):
            if len(reco_2dTracks['scifi'][p])!=3: raise ValueError('number of 2d tracks should be 3')
            for aTrack in reco_2dTracks['scifi'][p]:
                rho_scifi = {s*10+p:0 for s in range(1,6)} 
                rc = aTrack.Fit("pol1","QS")
                line = aTrack.GetFunction("pol1")
                for s in range(1,6):
                    mat  = 1
                    sipm = 1
                    channel = 64
                    plane_id=channel+1000*sipm+10000*mat+100000*p+1000000*s
                    self.geo.modules['Scifi'].GetSiPMPosition(plane_id,A,B)
                    ex = line.Eval((A[2]+B[2])/2)
                    for aHit in self.eventTree.Digi_ScifiHits :
                        if not aHit.isValid() : continue
                        detID = aHit.GetDetectorID()
                        so = detID//100000
                        if so !=s*10+p :continue
                        self.geo.modules['Scifi'].GetSiPMPosition(aHit.GetDetectorID(),A,B)
                        if p == 0 : pos=(A[1]+B[1])/2
                        else : pos=(A[0]+B[0])/2
                        if abs(pos-ex)<1. :
                            rho_scifi[s*10+p]+=1
                    if rho_scifi[s*10+p] > self.rho_cut:
                        noisy_stations[s]=1
        if sum(noisy_stations.values()) > 1 : passes = False
        return passes

class DSHitDensityCut(BaseCut):

    def __init__(self,rho):
        self.rho_cut=rho
        self.eventTree = tree
        self.geo = geo
 
    def pass_cut(self,reco_2dTracks):
        A,B = ROOT.TVector3(), ROOT.TVector3()
        noisy_stations = {s:0 for s in range(1,5)}
 
        sorted_hits = {s * 10 + o: [] for s in range(1, 5) for o in range(2)}
        for aHit in eventTree.Digi_MuFilterHits:
            if aHit.GetSystem() !=3 : continue
            detID = aHit.GetDetectorID()
            station  = (detID%10000)//1000  # station number
            bar = (detID%1000)
            so = (station+1)*10+(bar>59)
            sorted_hits[so].append(aHit)
        for p in range(2):
            for aTrack in reco_2dTracks['ds'][p]:
                rho_ds = {s*10+o:0 for s in range(1,5) for o in range(2)} 
                rc = aTrack.Fit("pol1","QS")
                line = aTrack.GetFunction("pol1")
                for so in sorted_hits:
                    if so == 40 : continue
                    if not so%10 == p : continue
                    geo.modules['MuFilter'].GetPosition(sorted_hits[so][0].GetDetectorID(),A,B)
                    ex = line.Eval((A[2]+B[2])/2)
                    for aHit in sorted_hits[so]:
                        geo.modules['MuFilter'].GetPosition(aHit.GetDetectorID(),A,B)
                        if p == 0 : pos=(A[1]+B[1])/2
                        else : pos=(A[0]+B[0])/2
                        if abs(pos-ex) < 3.1 : rho_ds[so]+=1
                    if rho_ds[s*10+p] > self.rho_cut:
                        noisy_stations[s]=1
                print(rho_ds)
        print(noisy_stations)
        if sum(noisy_stations.values()) > 1 : passes = False
        print(passes)
        return passes


class EventDeltaTCut(BaseCut):

    def __init__(self,delta_e, delta_t, eventTree):
        self.delta_t = delta_t
        
    def pass_cut(self):
        header=eventTree.EventHeader
        current_entry = self.eventTree.GetReadEvent()
        current_time  = header.GetEventTime()
        passes = True
        rc= self.eventTree.GetEvent(current_entry+delta_e)
        sign = (delta_e>0)-(delta_e<0)
        if -sign*(current_time-header.GetEventTime()) <= delta_t:
            passes = False
        plot_var = -sign*(current_time-header.GetEventTime())
        rc=self.eventTree.GetEvent(current_entry)
        return passes

class FiducialCut(BaseCut):

    def __init__(self):
        print("cut for 31x31 square cm fiducial region")
        
    def pass_cut(self,reco_2dTracks):
#        −42 cm ≤ x ≤ −11 cm and 18 cm ≤ y ≤ 49 cm
        passes = True
        for p in range(2):
            for aTrack in reco_2dTracks['scifi'][p]:
                rc = aTrack.Fit("pol1","QS")
                line = aTrack.GetFunction("pol1")
                pos=line.Eval(300.)
                if p==0:
                    if pos<18. or pos > 49:
                        passes =False
                else :
                    if pos<-42 or pos>-11:
                        passes = False
            if len(reco_2dTracks['scifi'][p])!=3:
                    raise ValueError('Number of 2d tracks are not 3 !')
        return passes







