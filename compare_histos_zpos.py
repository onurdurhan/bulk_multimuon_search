import ROOT

def init_style():

  #for the canvas:
  ROOT.gStyle.SetCanvasBorderMode(0)
  ROOT.gStyle.SetCanvasColor(ROOT.kWhite)
  ROOT.gStyle.SetCanvasDefH(600) #Height of canvas
  ROOT.gStyle.SetCanvasDefW(600) #Width of canvas
  ROOT.gStyle.SetCanvasDefX(10)   #Position on screen
  ROOT.gStyle.SetCanvasDefY(10)


  ROOT.gStyle.SetPadBorderMode(0)
  ROOT.gStyle.SetPadColor(ROOT.kWhite)
  ROOT.gStyle.SetPadGridX(False)
  ROOT.gStyle.SetPadGridY(False)
  ROOT.gStyle.SetGridColor(0)
  ROOT.gStyle.SetGridStyle(3)
  ROOT.gStyle.SetGridWidth(1)

  #For the frame:
  ROOT.gStyle.SetFrameBorderMode(0)
  ROOT.gStyle.SetFrameBorderSize(1)
  ROOT.gStyle.SetFrameFillColor(0)
  ROOT.gStyle.SetFrameFillStyle(0)
  ROOT.gStyle.SetFrameLineColor(1)
  ROOT.gStyle.SetFrameLineStyle(1)
  ROOT.gStyle.SetFrameLineWidth(3)

  #For the histo:
  ROOT.gStyle.SetHistLineColor(ROOT.kBlue)
  ROOT.gStyle.SetHistLineStyle(0)
  ROOT.gStyle.SetHistLineWidth(2)


  ROOT.gStyle.SetEndErrorSize(2)


  ROOT.gStyle.SetMarkerStyle(20)

  #For the fit/function:
  ROOT.gStyle.SetOptFit(1)
  ROOT.gStyle.SetFitFormat("5.4g")
  ROOT.gStyle.SetFuncColor(2)
  ROOT.gStyle.SetFuncStyle(1)
  ROOT.gStyle.SetFuncWidth(1)

  #For the date:
  ROOT.gStyle.SetOptDate(0)


  # For the statistics box:
  ROOT.gStyle.SetOptFile(0)
  ROOT.gStyle.SetOptStat(0) # To display the mean and RMS:   SetOptStat("mr")
  ROOT.gStyle.SetStatColor(ROOT.kWhite)
  ROOT.gStyle.SetStatFont(42)
  ROOT.gStyle.SetStatFontSize(0.025)
  ROOT.gStyle.SetStatTextColor(1)
  ROOT.gStyle.SetStatFormat("6.4g")
  ROOT.gStyle.SetStatBorderSize(1)
  ROOT.gStyle.SetStatH(0.1)
  ROOT.gStyle.SetStatW(0.15)

  # Margins:
  ROOT.gStyle.SetPadTopMargin(0.05)
  ROOT.gStyle.SetPadBottomMargin(0.15)
  ROOT.gStyle.SetPadLeftMargin(0.15)
  ROOT.gStyle.SetPadRightMargin(0.05)

  # For the Global title:

#  ROOT.gStyle.SetOptTitle(0)
  ROOT.gStyle.SetTitleFont(42)
  ROOT.gStyle.SetTitleColor(1)
  ROOT.gStyle.SetTitleTextColor(1)
  ROOT.gStyle.SetTitleFillColor(10)
  ROOT.gStyle.SetTitleFontSize(0.05)

  # For the axis titles:

  ROOT.gStyle.SetTitleColor(1, "XYZ")
  ROOT.gStyle.SetTitleFont(42, "XYZ")
  ROOT.gStyle.SetTitleSize(0.05, "XYZ")
  ROOT.gStyle.SetTitleXOffset(1)
  ROOT.gStyle.SetTitleYOffset(1.25)

  # For the axis labels:
    
  ROOT.gStyle.SetLabelColor(1, "XYZ")
  ROOT.gStyle.SetLabelFont(42, "XYZ")
  ROOT.gStyle.SetLabelOffset(0.007, "XYZ")
  ROOT.gStyle.SetLabelSize(0.04, "XYZ")
  

  # For the axis:

  ROOT.gStyle.SetAxisColor(1, "XYZ")
  ROOT.gStyle.SetStripDecimals(True)
  ROOT.gStyle.SetTickLength(0.03, "XYZ")
  ROOT.gStyle.SetNdivisions(510, "XYZ")
  ROOT.gStyle.SetPadTickX(1)  # To get tick marks on the opposite side of the frame
  ROOT.gStyle.SetPadTickY(1)
  

  # Change for log plots:
  #ROOT.gStyle.SetOptLogx(0)
  #ROOT.gStyle.SetOptLogy(0)
  #ROOT.gStyle.SetOptLogz(0)

  #Legend options:
  ROOT.gStyle.SetLegendBorderSize(0)
  ROOT.gStyle.SetLegendTextSize(0.022)

  # Postscript options:
  #ROOT.gStyle.SetPaperSize(20.,26.)
  #ROOT.gStyle.SetHatchesLineWidth(5)
  #ROOT.gStyle.SetHatchesSpacing(0.05)

def writeSND(pad,
  text_factor=0.9,
  text_offset=0.01,
  extratext="",
  text_in=False,
  rfrac=0,
  maintext='SND@LHC'):

  pad.Update()

  l = pad.GetLeftMargin()
  t = pad.GetTopMargin()
  r = pad.GetRightMargin()
  b = pad.GetBottomMargin()

  SNDTextSize = t*text_factor
  SNDTextVerticalOffset = text_offset

  pad.cd()

  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextAngle(0)
  latex.SetTextColor(ROOT.kBlack)

  latex.SetTextFont(61)
  latex.SetTextAlign(11)
  latex.SetTextSize(SNDTextSize)
  latex.SetText(0,0,maintext)

  sndX = SNDTextSize*2*(1-rfrac)

  if not text_in: latex.DrawLatex(l, 1-t+SNDTextVerticalOffset, maintext)
  else: latex.DrawLatex(l+0.03, 1-t-SNDTextVerticalOffset-1.2*SNDTextSize, maintext)

  extraTextSize = SNDTextSize*0.8
  latex.SetTextFont(52)
  latex.SetTextSize(extraTextSize)
  latex.SetTextAlign(11)
  if not text_in: latex.DrawLatex(l+0.03 + 1.5*sndX, 1-t+SNDTextVerticalOffset, extratext)
  else: latex.DrawLatex(l+0.03, 1-t-SNDTextVerticalOffset-2*SNDTextSize, extratext)
  pad.Update()
  return



rc=init_style()
# Load ROOT files
data_file = ROOT.TFile.Open("histos_data.root")  # Replace with actual data file
mc_file = ROOT.TFile.Open("histos_MC_Molasse.root")   #("histos_MC_bathces2.root")  # Replace with actual MC file
# Known luminosities
L_data = 4.827999999999999 #0.025  # Replace with actual data luminosity in fb^-1
L_MC   = 1. # Given MC luminosity in fb^-1
scale_factor = L_data / L_MC

# Histogram names
hist_data = [
    data_file.Get("data_histo_dmax_vs_dmin_0"),
    data_file.Get("data_histo_dmax_vs_dmin_1"),
]

# Define styles for different vertex position ranges
position_colors = {
    "7150_5000": ROOT.kRed,  # Solid line, width 2
    "5000_3000": ROOT.kBlue,  # Dashed line, width 2
    "3000_1000": ROOT.kGreen,  # Dotted line, width 2
    "1000_0"   : ROOT.kMagenta,  # Dash-dot line, width 2
}

good_colors = [
    ROOT.kBlue+1,      # Deep Blue
    ROOT.kRed-4,       # Dark Red
    ROOT.kGreen+2,     # Bright Green
    ROOT.kOrange+7,    # Strong Orange
    ROOT.kViolet-6,    # Soft Purple
    ROOT.kCyan+2,      # Bright Cyan
    ROOT.kMagenta-7,   # Deep Magenta
    ROOT.kYellow+1,    # Warm Yellow
    ROOT.kAzure+7,     # Sky Blue
    ROOT.kPink-3,      # Soft Pink
    ROOT.kSpring-9,    # Light Green
    ROOT.kTeal-6       # Muted Teal
]
#hs={ 0:ROOT.THStack("hs0"," dmax vs dmin"),1:ROOT.THStack("hs1"," dmax vs dmin")}
# MC histogram dictionary
RATIO=False
if RATIO:
    hist_data=[data_file.Get("data_dmin_dmax_ratio_0"),data_file.Get("data_dmin_dmax_ratio_1")]
    colors={}
    i=0
    for prod in ["MuonToMuonPair","GammaToMuonPair","PositronToMuonPair"]:
        colors[prod]={}
        for pos in ["7150_5000", "5000_3000","3000_1000","1000_0"]:
            colors[prod][pos]=good_colors[i]
            i+=1
else :
    colors={}
    i=0
    for prod in ["MuonToMuonPair","GammaToMuonPair","PositronToMuonPair"]:
        colors[prod]={}
        for pos in ["7150_5000", "5000_3000","3000_1000","1000_0"]:
            colors[prod][pos]=good_colors[i]
            i+=1




mc_hists = {}
hs={}
signal_bkg = {}
for view in range(2):
    hs[view]={"MuonToMuonPair":{},"GammaToMuonPair":{},"PositronToMuonPair":{}}
    mc_hists[view] = {"MuonToMuonPair":{}, "GammaToMuonPair":{}, "PositronToMuonPair":{}}
    signal_bkg[view] = {'Signal':{},'Background':{}}
    if RATIO:
        hs[view] = ROOT.THStack(f"hs_{view}"," dmin dmax ratio")
    cloned={"7150_5000":False, "5000_3000":False, "3000_1000":False, "1000_0":False}

    for prod in mc_hists[view] :
        mc_hists[view][prod]={"7150_5000":{}, "5000_3000":{}, "3000_1000":{}, "1000_0":{}}
        if prod=='MuonToMuonPair':
            signal_bkg[view]['Signal'] = mc_hists[view][prod]
        else:
            if prod not in mc_hists[view].keys() : signal_bkg[view]['Background'] = mc_hists[view][prod]
        for pos_range in mc_hists[view][prod]:
            if RATIO:
                hist_name = f"{prod}_dmin_dmax_ratio_{view}_{pos_range}"
            else:
                hist_name = f"{prod}_histo_dmax_vs_dmin_{view}_{pos_range}"
            if prod == "MuonToMuonPair":
                mc_hists[view][prod][pos_range] = mc_file.Get(hist_name)
                signal_bkg[view]['Signal'][pos_range] = mc_hists[view][prod][pos_range]
            else :
                if cloned[pos_range] :
                    print("trying to add", hist_name, " histo on", signal_bkg[view]['Background'][pos_range], "for ",pos_range)
                    signal_bkg[view]['Background'][pos_range].Add(mc_file.Get(hist_name))
                    print('Added histograms for bckg', mc_file.Get(hist_name))
                else :
                    signal_bkg[view]['Background'][pos_range] = mc_file.Get(hist_name).Clone(f'{view}_{prod}_{pos_range}')
                    print('cloned histos ', signal_bkg[view]['Background'][pos_range])
                    cloned[pos_range] = True
    for prod in ['Signal','Background']:
        for pos_range in mc_hists[view]['MuonToMuonPair'].keys():
            print(signal_bkg[view][prod][pos_range])




# Define dmax slice bins
num_slices = 3 # Adjust based on data range
canvases = [ROOT.TCanvas(f"c_{idx}", f"Comparison of Data and MC {idx}",1200,800) for idx in range(2)]
for canvas in canvases:
    if not RATIO:
        canvas.Divide(2, 2)  # Divide canvas into pads

legends=[]
for view, canvas in enumerate(canvases):
    if RATIO:
        pad=canvas.cd()
        legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
        max_data = hist_data[view].GetMaximum() + hist_data[view].GetBinError(hist_data[view].GetMaximumBin())
        for prod in signal_bkg[view]:
#            if prod == 'Signal' :continue
            for pos_range in signal_bkg[view][prod]:
                hist=signal_bkg[view][prod][pos_range]
                hist=hist.Rebin(3)
                max_mc = hist.GetMaximum() + hist.GetBinError(hist.GetMaximumBin())
                maxY=max(max_data, max_mc)
#                hist.SetLineColor(colors[prod][pos_range])
#                hist.SetFillColor(colors[prod][pos_range])
#                hist.SetFillStyle(1001)
#                hist.SetFillColorAlpha(colors[prod][pos_range],0.45)

                if prod=='Signal':
                    hist.SetFillColor(colors['MuonToMuonPair'][pos_range])
                    hist.SetLineColor(colors['MuonToMuonPair'][pos_range])
                else:
                    hist.SetFillColor(colors['GammaToMuonPair'][pos_range])
                    hist.SetLineColor(colors['GammaToMuonPair'][pos_range])
                hist.SetLineWidth(2)
                hist.Scale(scale_factor)  # Apply luminosity scaling and normalize
                hist.GetXaxis().SetTitle("#frac{dmin}{dmax} [cm]")
                hs[view].Add(hist)
                legend.AddEntry(hist, prod+pos_range, "f")
            rc=hs[view].SetTitle("dmin dmax ratio ")
#        hs[view].SetMaximum(3 * maxY)
        hs[view].Draw("HIST")
        hs[view].GetXaxis().SetTitle("#frac{d_{min}}{d_{max}}")
        hist_data[view]=hist_data[view].Rebin(3)
        hist_data[view].SetLineColor(ROOT.kBlack)
        hist_data[view].SetMarkerStyle(20)
        hist_data[view].Draw("SAME E")
        legend.AddEntry(hist_data[view], "Data", "lep")
        legend.Draw("SAME")
        legends.append(legend)
        rc=writeSND(pad)
        canvas.SaveAs(f"comparison_slices_{view}.png")
        continue

    else :
        print("switch to slices")
    for i in range(num_slices):
        pad = canvas.cd(i + 1)
        legend = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
        if i < 2 :
            n_first =  i
            n_last = i+1
        else :
            n_first = i
            n_last = hist_data[view].GetNbinsX()
        dmin_min = hist_data[view].GetXaxis().GetBinLowEdge(i+1)
        dmin_max = hist_data[view].GetXaxis().GetBinUpEdge(n_last) 
        proj_data = hist_data[view].ProjectionY(f"data_proj_{view}_{i}", n_first,n_last)
        proj_data.Rebin(1)
        max_data = proj_data.GetMaximum() + proj_data.GetBinError(proj_data.GetMaximumBin());
        hs[view][n_last]=ROOT.THStack(f"hs_{view}_{n_last}"," dmax vs dmin")
        for prod  in signal_bkg[view]:#mc_hists[view]:
#            if prod =='Signal':continue
            for pos_range in signal_bkg[view][prod]:#mc_hists[view][prod]:
                hist=signal_bkg[view][prod][pos_range]#mc_hists[view][prod][pos_range]
                proj_mc =hist.ProjectionY(f"mc_{prod}_{view}_{i}_{pos_range}", n_first, n_last)
                proj_mc.Rebin(1)
                max_mc = proj_mc.GetMaximum() + proj_mc.GetBinError(proj_mc.GetMaximumBin())
                maxY=max(max_data, max_mc)
#                proj_mc.SetLineColor(position_colors[pos_range])
                if prod=='Signal':
                    proj_mc.SetFillColor(colors['MuonToMuonPair'][pos_range])
                    proj_mc.SetLineColor(colors['MuonToMuonPair'][pos_range])
                else:
                    proj_mc.SetFillColor(colors['GammaToMuonPair'][pos_range])
                    proj_mc.SetLineColor(colors['GammaToMuonPair'][pos_range])
                proj_mc.SetLineWidth(2)
                proj_mc.Scale(scale_factor)  # Apply luminosity scaling and normalize
                hs[view][n_last].Add(proj_mc)
                legend.AddEntry(proj_mc, prod+pos_range, "f")
        rc=hs[view][n_last].SetTitle(f"{dmin_min} < dmin < {dmin_max}")
        hs[view][n_last].SetMaximum(1.5 * maxY)
        hs[view][n_last].Draw("HIST")
        hs[view][n_last].GetXaxis().SetTitle("dmax [cm]")
        proj_data.SetLineColor(ROOT.kBlack)
        proj_data.SetMarkerStyle(20)
        proj_data.Draw("SAME E")
        legend.AddEntry(proj_data, "Data", "lep")
        legend.Draw("SAME")
        legends.append(legend)
        rc=writeSND(pad)

#    canvas.Update()
    canvas.SaveAs(f"comparison_slices_{view}.pdf")
