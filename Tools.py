import ROOT
class Tools:
    def __init__(self,trackTask,geo):
        self.trackTask=trackTask
        self.geo = geo
        print("auxillary tools initiated")
        
    def scifi_clusters(self):
        self.trackTask.clusScifi.Clear()
        self.trackTask.scifiCluster()
        clusters = self.trackTask.clusScifi
        return clusters

    def scifi_clusters_per_plane(self,clusters):
        sortedClusters = {s * 10 + o: [] for s in range(1, 6) for o in range(2)}
        for aCl in clusters:
            so = aCl.GetFirst()//100000
            sortedClusters[so].append(aCl)
        return sortedClusters
##### get average plane pos


    def getAverageZpositions(self):
        A,B = ROOT.TVector3(), ROOT.TVector3()
        zPos={'MuFilter':{},'Scifi':{}}
        Scifi  = self.geo.modules['Scifi']
        MuFilter  = self.geo.modules['MuFilter']
        nScifi = Scifi.GetConfParI('Scifi/nscifi')
        for s in self.systemAndPlanes():
            for plane in range(self.systemAndPlanes()[s]):
                bar = 4
                p = plane
                if s==3 and (plane%2==0 or plane==7):
                    bar = 90
                    p = plane//2
                elif s==3 and plane%2==1:
                    bar = 30
                    p = plane//2
                MuFilter.GetPosition(s*10000+p*1000+bar,A,B)
                zPos['MuFilter'][s*10+plane] = (A.Z()+B.Z())/2.
        for s in range(1,nScifi+1):
            mat  = 1
            sipm = 1
            channel = 64
            for o in range(2):
                Scifi.GetPosition(channel+1000*sipm+10000*mat+100000*o+1000000*s,A,B)
                zPos['Scifi'][s*10+o] = (A.Z()+B.Z())/2.
        return zPos


    def systemAndOrientation(self,s,plane):
        if s==1 or s==2: return "horizontal"
        if plane%2==1 or plane == 6: return "vertical"
        return "horizontal"

    def systemAndPlanes(self):
        MuFilter = self.geo.modules['MuFilter']
        systemAndPlanes = {1:MuFilter.GetConfParI("MuFilter/NVetoPlanes"),
                           2:MuFilter.GetConfParI("MuFilter/NUpstreamPlanes"),
                           3:2*MuFilter.GetConfParI("MuFilter/NDownstreamPlanes")-1} # to arrive at 7 DS planes
        return systemAndPlanes
    def systemAndBars(self) :   
        MuFilter = self.geo.modules['MuFilter']
        systemAndBars =   {1:MuFilter.GetConfParI("MuFilter/NVetoBars"),
                           2:MuFilter.GetConfParI("MuFilter/NUpstreamBars"),
                           3:MuFilter.GetConfParI("MuFilter/NDownstreamBars")}
        return systemAndBars
