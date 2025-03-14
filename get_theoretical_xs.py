#!/usr/bin/env python
import os
import sys
import ROOT
import numpy as np
import shipunit as u
import shipRoot_conf
import rootUtils as ut
from ShipGeoConfig import ConfigRegistry
from argparse import ArgumentParser
# Load required GEANT4 libraries
ROOT.gSystem.Load("libG4processes")
ROOT.gSystem.Load("libG4tracking")
ROOT.gSystem.Load("libG4track")
ROOT.gSystem.Load("libG4materials")
ROOT.gSystem.Load("libG4particles")
ROOT.gSystem.Load("libG4global")
ROOT.gSystem.Load("libG4physicslists")
ROOT.gSystem.Load("libG4parmodels")
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuonToMuonPairProductionModel.hh"')




mcEngine     = "TGeant4"

inactivateMuonProcesses  = False

MCTracksWithHitsOnly                  = False  # copy particles which produced a hit and their history
MCTracksWithEnergyCutOnly    = True # copy particles above a certain kin energy cut
MCTracksWithHitsOrEnergyCut = False # or of above, factor 2 file size increase compared to MCTracksWithEnergyCutOnly

parser = ArgumentParser()

parser.add_argument("--Ntuple",  dest="ntuple",  help="Use ntuple as input", required=False, action="store_true")
parser.add_argument('--eMin', type=float, help="energy cut", dest='ecut', default=-1.)
parser.add_argument('--zMax', type=float, help="max distance to apply energy cut", dest='zmax', default=70000.)
parser.add_argument("-n", "--nEvents",dest="nEvents",  help="Number of events to generate", required=False,  default=100, type=int)
parser.add_argument("-i", "--firstEvent",dest="firstEvent",  help="First event of input file to use", required=False,  default=0, type=int)
parser.add_argument("-s", "--seed",dest="theSeed",  help="Seed for random number. Only for experts, see TRrandom::SetSeed documentation", required=False,  default=0, type=int)
parser.add_argument("-f", dest="inputFile", help="Input file if not default file", required=False, default="LHC_-160urad_magfield_2022TCL6_muons_rock_2e8pr_30GeVcut.root")
parser.add_argument("-o", "--output",dest="outputDir",  help="Output directory", required=False,  default=".")
parser.add_argument("--boostFactor", dest="boostFactor",  help="boost mu brems", required=False, type=float,default=0)
parser.add_argument("--debug",   dest="debug",   help="debugging mode, check for overlaps", required=False, action="store_true")
parser.add_argument("-D", "--display", dest="eventDisplay", help="store trajectories", required=False, action="store_true")
parser.add_argument("--EmuDet","--nuTargetActive",dest="nuTargetPassive",help="activate emulsiondetector", required=False,action="store_false")
parser.add_argument("--NagoyaEmu","--useNagoyaEmulsions",dest="useNagoyaEmulsions",help="use bricks of 57 Nagoya emulsion films instead of 60 Slavich", required=False,action="store_true")
parser.add_argument("-y", dest="year", help="specify the year to generate the respective TI18 detector setup", required=False, type=int, default=2024)

options = parser.parse_args()

simEngine = "Ntuple"

inputFile = options.inputFile

hist={}
ut.bookHist(hist,"xs","",1000,0,300,1000,0,1500.)






print("SND@LHC setup for",simEngine,"to produce",options.nEvents,"events")

ROOT.gRandom.SetSeed(options.theSeed)  # this should be propagated via ROOT to Pythia8 and Geant4VMC
shipRoot_conf.configure(0)     # load basic libraries, prepare atexit for python

snd_geo = ConfigRegistry.loadpy("$SNDSW_ROOT/geometry/sndLHC_geom_config.py",
                                   nuTargetPassive = options.nuTargetPassive,
                                   useNagoyaEmulsions = options.useNagoyaEmulsions,
                                   year=options.year)

tag = simEngine+"-"+mcEngine

if not os.path.exists(options.outputDir):
  os.makedirs(options.outputDir)
tag+='_muToMuPair'
outFile = "%s/sndLHC.%s.root" % (options.outputDir, tag)

# rm older files !!! 
for x in os.listdir(options.outputDir):
  if not x.find(tag)<0: os.system("rm %s/%s" % (options.outputDir, x) )
# Parameter file name
parFile="%s/ship.params.%s.root" % (options.outputDir, tag)

# In general, the following parts need not be touched, except for user task
# ========================================================================

# -----Timer--------------------------------------------------------
timer = ROOT.TStopwatch()
timer.Start()
# ------------------------------------------------------------------------
# -----Create simulation run----------------------------------------
run = ROOT.FairRunSim()
run.SetName(mcEngine)  # Transport engine
run.SetSink(ROOT.FairRootFileSink(outFile))  # Output file
run.SetUserConfig("g4Config.C") # user configuration file default g4Config.C 
rtdb = run.GetRuntimeDb() 
# -----Create geometry----------------------------------------------
import shipLHC_conf as sndDet_conf
modules = sndDet_conf.configure(run,snd_geo)

# -----Create PrimaryGenerator--------------------------------------
primGen = ROOT.FairPrimaryGenerator()

ut.checkFileExists(inputFile)
Ntuplegen = ROOT.NtupleGenerator_FLUKA()
Ntuplegen.SetZ(snd_geo.Floor.z)
Ntuplegen.Init(inputFile,options.firstEvent)
primGen.AddGenerator(Ntuplegen)
#   options.nEvents = Ntuplegen.GetNevents()


if options.ecut > 0:   
   modules['Floor'].SetEmin(options.ecut)
   modules['Floor'].SetZmax(options.zmax)

#
run.SetGenerator(primGen)
# ------------------------------------------------------------------------

#---Store the visualiztion info of the tracks, this make the output file very large!!
#--- Use it only to display but not for production!
if options.eventDisplay: run.SetStoreTraj(ROOT.kTRUE)
else:            run.SetStoreTraj(ROOT.kFALSE)



# -----Initialize simulation run------------------------------------
run.Init()
gMC = ROOT.TVirtualMC.GetMC()
fStack = gMC.GetStack()
if MCTracksWithHitsOnly:
 fStack.SetMinPoints(1)
 fStack.SetEnergyCut(-100.*u.MeV)
elif MCTracksWithEnergyCutOnly:
 fStack.SetMinPoints(-1)
 fStack.SetEnergyCut(100.*u.MeV)
elif MCTracksWithHitsOrEnergyCut: 
 fStack.SetMinPoints(1)
 fStack.SetEnergyCut(100.*u.MeV)
elif options.deepCopy: 
 fStack.SetMinPoints(0)
 fStack.SetEnergyCut(0.*u.MeV)


ROOT.gROOT.ProcessLine('#include "Geant4/G4ProcessTable.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuBremsstrahlung.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4GammaConversionToMuons.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuPairProduction.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4AnnihiToMuPair.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuonToMuonPairProduction.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuonPlus.hh"')
ROOT.gROOT.ProcessLine('#include "Geant4/G4MuonMinus.hh"')

gProcessTable = ROOT.G4ProcessTable.GetProcessTable()
 # only muon interaction
 # procBrems        = gProcessTable.FindProcess(ROOT.G4String('muBrems'),ROOT.G4String('mu+'))
 # muPairProd = gProcessTable.FindProcess(ROOT.G4String('muPairProd'),ROOT.G4String('mu+'))
 # muPairProd.SetCrossSectionBiasingFactor(options.boostFactor)
 # procBrems.SetCrossSectionBiasingFactor(options.boostFactor)
 # muon pair production
gammaToMuPair = gProcessTable.FindProcess(ROOT.G4String('GammaToMuPair'),ROOT.G4String('gamma'))
#gammaToMuPair.SetCrossSecFactor(options.boostFactor) 
AnnihiToMuPair = gProcessTable.FindProcess(ROOT.G4String('AnnihiToMuPair'),ROOT.G4String('e+'))
#AnnihiToMuPair.SetCrossSecFactor(options.boostFactor)
MuonToMuonPair = gProcessTable.FindProcess(ROOT.G4String('muToMuonPairProd'),ROOT.G4String('mu+'))
#MuonToMuonPair.SetCrossSectionBiasingFactor(100)

mygMC = ROOT.TGeant4.GetMC()
# -----Start run----------------------------------------------------
run.Run(1)
# -----Runtime database---------------------------------------------
"""
kParameterMerged = ROOT.kTRUE
parOut = ROOT.FairParRootFileIo(kParameterMerged)
parOut.open(parFile)
rtdb.setOutput(parOut)
rtdb.saveOutput()
rtdb.printParamContexts()
getattr(rtdb,"print")()
# ------------------------------------------------------------------------
geoFile = "%s/geofile_full.%s.root" % (options.outputDir, tag)
run.CreateGeometryFile(geoFile)
# save detector parameters dictionary in geofile
import saveBasicParameters
saveBasicParameters.execute(geoFile,snd_geo)

# ------------------------------------------------------------------------
"""

fluka_file=ROOT.TFile(inputFile)
sTree = fluka_file.nt
#nist_manager=ROOT.G4NistManager.Instance()
#muon = ROOT.G4MuonPlus.MuonPlusDefinition()
muon = ROOT.G4MuonMinus.MuonMinusDefinition()
rock = ROOT.G4Material.GetMaterial("MolasseRock")
tungsten = ROOT.G4Material.GetMaterial("tungstensifon")
print("rock number density ", rock.GetTotNbOfAtomsPerVolume())
#Pb=nist_manager.FindOrBuildElement("Pb")
#lead = ROOT.G4Material("lead",11.34, 1)
#lead.AddElement(Pb, 1)
#couple = ROOT.G4MaterialCutsCouple(material)
model=MuonToMuonPair.EmModel()
mm2_to_nb = 1E31
array = np.arange(0, 10001)


# Arrays to accumulate cross sections and weights
cross_sections = []
weights = []

# Loop over the tree and calculate the cross section for each energy value
for event in sTree:
    energy = sTree.E*1E3
    if energy > 300*1E3 : continue
    # Calculate the cross section for the current energy
    cross_section = (model.CrossSectionPerVolume(tungsten, muon, energy)/rock.GetTotNbOfAtomsPerVolume())*mm2_to_nb
    # Accumulate the cross sections and weights (assuming weight = 1 for raw counts)
    cross_sections.append(cross_section)
    weights.append(sTree.w)  # Use event weight if available in the tree
    hist["xs"].Fill(sTree.E, cross_section)


# Calculate the weighted average cross section

if np.sum(weights) > 0:
    average_cross_section = np.sum(np.array(cross_sections) * np.array(weights)) / np.sum(weights)
    average_cross_section_nb = average_cross_section
    print(f"Weighted Average Cross Section: {average_cross_section_nb:.6f} nb")
else:
    print("No valid energy data found in the tree!")









# -----Finish-------------------------------------------------------
timer.Stop()
rtime = timer.RealTime()
ctime = timer.CpuTime()
print(' ') 
print("Macro finished succesfully.") 

print("Output file is ",  outFile) 
#print("Geometry file is ",geoFile)
print("Real time ",rtime, " s, CPU time ",ctime,"s")

# ------------------------------------------------------------------------
def checkOverlaps(removeScifi=False):
 sGeo = ROOT.gGeoManager
 if removeScifi:
    for n in range(1,6):
       Hscifi = sGeo.FindVolumeFast('ScifiVolume'+str(n))
       removalList = []
       for x in Hscifi.GetNodes():
             if x.GetName().find('Scifi')==0: removalList.append(x)
       for x in removalList: Hscifi.RemoveNode(x)
 sGeo.SetNmeshPoints(10000)
 sGeo.CheckOverlaps(0.1)  # 1 micron takes 5minutes
 sGeo.PrintOverlaps()
# check subsystems in more detail
 for x in sGeo.GetTopNode().GetNodes(): 
   x.CheckOverlaps(0.0001)
   sGeo.PrintOverlaps()

def checkOverlapsWithGeant4():
 # after /run/initialize, but prints warning messages, problems with TGeo volume
 mygMC = ROOT.TGeant4.GetMC()
 mygMC.ProcessGeantCommand("/geometry/test/recursion_start 0")
 mygMC.ProcessGeantCommand("/geometry/test/recursion_depth 2")
 mygMC.ProcessGeantCommand("/geometry/test/run")
