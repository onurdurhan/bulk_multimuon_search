import subprocess
from argparse import ArgumentParser
from MuonTridentEngine import  MuonTridentEngine
# Setup argument parser
parser = ArgumentParser()
parser.add_argument("operation", choices=["filter", "search"], help="operation to perform: filter or search")
parser.add_argument("-f", "--inputFile", dest="inputFile", help="input file", required=True)
parser.add_argument("-g", "--geoFile", dest="geoFile", help="geo file", required=False)
parser.add_argument("-o", "--outFile", dest="outFile", help="output file", required=True)
parser.add_argument("-r", "--runNumber", dest="runNumber", help="run number", required=False)
parser.add_argument("-p", "--path", dest="path",help="path to physics data", required=False)
options = parser.parse_args()
if not options.geoFile:
    if options.path.find('2022')>0 :
        if int(options.runNumber) < 4575:
            options.geoFile =  options.path+"/geofile_sndlhc_TI18_V3_08August2022.root"
        elif int(options.runNumber) <= 4855:
            options.geoFile =  options.path+"/geofile_sndlhc_TI18_V5_14August2022.root"
        elif int(options.runNumber) < 5172:
            options.geoFile =  options.path+"/geofile_sndlhc_TI18_V6_08October2022.root"
        else:
            options.geoFile =  options.path+"/geofile_sndlhc_TI18_V7_22November2022.root"
    elif options.path.find('2023')>0 : 
        options.geoFile =  options.path+"/geofile_sndlhc_TI18_V4_2023.root"
    else:
        print("No such run is taken yet")
        sys.exit()








# Function to run the filtering operation
def run_filter(options):
#    parts = options.inputFile.split('/')
#    fname = parts[-1]
#    partition = int(fname.split('-')[-1].split('.')[0])
#    options.outFile = "filtered_sndsw_run_"+str(options.runNumber)+"_"+str(partition)+".root"
    engine = MuonTridentEngine(options)
    rc = engine.filter_time()

# Function to run the search operation
def run_search(options):
#    parts = options.inputFile.split('/')
#    fname = options.inputFile #parts[-1]
#    partition = int(fname.split('_')[-1].split('.')[0])
#    options.inputFile =  "filtered_sndsw_run_"+str(options.runNumber)+"_"+str(partition)+".root"
#    options.outFile = "mu3_sndsw_run_"+str(options.runNumber)+"_"+str(partition)+".root"
    engine = MuonTridentEngine(options)
    rc = engine.search_mu3()

# Main execution logic based on operation
if options.operation == "filter":
    run_filter(options)
elif options.operation == "search":
    run_search(options)

