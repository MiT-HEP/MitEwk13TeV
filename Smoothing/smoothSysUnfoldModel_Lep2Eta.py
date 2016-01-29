import ROOT
import GausKernelSmoother
import LatexDocument
import math

def SetPlotStyle():
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat()
    ROOT.gStyle.SetOptFit(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetFrameLineWidth(1)
    ROOT.gStyle.SetPadBottomMargin(0.16)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadRightMargin(0.05)

    ROOT.gStyle.SetLabelSize(0.05,"X")
    ROOT.gStyle.SetLabelSize(0.05,"Y")
    ROOT.gStyle.SetLabelOffset(0.01,"Y")
    ROOT.gStyle.SetTickLength(0.04,"X")
    ROOT.gStyle.SetTickLength(0.04,"Y")
    ROOT.gStyle.SetLineWidth(1)
    ROOT.gStyle.SetTickLength(0.04 ,"Z")

    ROOT.gStyle.SetTitleSize(0.1)
    ROOT.gStyle.SetTitleSize(0.05,"X")
    ROOT.gStyle.SetTitleSize(0.05,"Y")
    ROOT.gStyle.SetTitleOffset(1.4,"X")
    ROOT.gStyle.SetTitleOffset(1.7,"Y")
    ROOT.gStyle.SetTitleFont(42, "X")
    ROOT.gStyle.SetTitleFont(42, "Y")
    ROOT.gStyle.SetLabelFont(42, "X")
    ROOT.gStyle.SetLabelFont(42, "Y")

    ROOT.gStyle.SetLegendBorderSize(0); 
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(1)
    ROOT.gROOT.ForceStyle()


def copyHistos(fileIn,fileOut):
    listOfDir = []
    listOfHisto = []
    next = ROOT.TIter(fileIn.GetListOfKeys())
    key = next()
    while key:
        if key.GetClassName().find("TDirectory")==0:
            listOfDir.append(key.GetName())
        elif key.GetClassName()=="TH1D":
            listOfHisto.append(key.GetName())
        key = next()

    fileOut.cd()
    for histName in listOfHisto:
        histo = fileIn.Get(histName)
        histo.__class__ = ROOT.TH1D
        histoCopy = histo.Clone()
        histoCopy.__class__ = ROOT.TH1D
        histoCopy.SetDirectory(0)
        histoCopy.Write()

    for dir in listOfDir:
        fileOut.mkdir(dir)
        fileOut.cd(dir)
        for histName in listOfHisto:
            histo = fileIn.Get(dir+"/"+histName)
            if not histo:
                raise StandardError("ERROR: cannot find "+dir+"/"+histName)
            histo.__class__ = ROOT.TH1D
            histoCopy = histo.Clone()
            histoCopy.__class__ = ROOT.TH1D
            histoCopy.SetDirectory(0)
            histoCopy.Write()

def createInclusiveCrossSection(exclHisto):
    name = exclHisto.GetName().replace("Excl", "Incl")
    inclHisto = exclHisto.Clone(name)
    inclHisto.__class__ = ROOT.TH1D
    for b in range(1,inclHisto.GetNbinsX()+1):
        sum = 0.
        error = 0.
        for b2 in range(b,exclHisto.GetNbinsX()+2):
            sum += exclHisto.GetBinContent(b2)
            error += exclHisto.GetBinError(b2)*exclHisto.GetBinError(b2)
        inclHisto.SetBinContent(b, sum)
        inclHisto.SetBinError(b,math.sqrt(error))
    return inclHisto

def getBinomialError(b1,b2,e1,e2):
    # binomial error of b1/b2
    # see http://root.cern.ch/root/html532/src/TH1.cxx.html#hh46dE 
    error = 0.
    if b1!=b2:
        w = b1/b2
        error = math.sqrt(abs( ( (1.-2.*w)*e1*e1 + w*w*e2*e2 )/(b2*b2) ))
    return error

def createMultiplicityRatio(multHisto):
    name = multHisto.GetName().replace("Excl","ExclRatio").replace("Incl","InclRatio")
    ratioHisto = multHisto.Clone(name)
    ratioHisto.__class__ = ROOT.TH1D
    nbins = multHisto.GetNbinsX()
    for b in range(2,nbins+1):
        nPlus1 = multHisto.GetBinContent(b)
        n = multHisto.GetBinContent(b-1)
        nPlus1E = multHisto.GetBinError(b)
        nE = multHisto.GetBinError(b-1)
        ratio = 0.
        ratioE = 0.
        if n>0.:
            ratio = nPlus1/n
            if "Incl" in name:
                # get binomial error
                ratioE = getBinomialError(nPlus1,n,nPlus1E,nE)
            else:
                ratioE = math.sqrt((nPlus1*nE)**2. + (n*nPlus1E)**2.)/(n*n)
        ratioHisto.SetBinContent(b, ratio)
        ratioHisto.SetBinError(b, ratioE)
    ratioHisto.SetBinContent(1, 0.)
    ratioHisto.SetBinError(1, 0.)
    return ratioHisto

def createInclusiveAndRatio(file, directory):
    print "Creating inclusive multiplicities and ratios in", directory
    listOfHisto = []
    next = ROOT.TIter(file.Get(directory).GetListOfKeys())
    key = next()
    while key:
        if key.GetClassName()=="TH1D" and "nExcl" in key.GetName():
            listOfHisto.append(key.GetName())
        key = next()
    file.cd(directory)
    for histName in listOfHisto:
        hist = file.Get(directory+"/"+histName)
        hist.__class__ = ROOT.TH1D
        inclHist = createInclusiveCrossSection(hist)
        exclRatio = createMultiplicityRatio(hist)
        inclRatio =createMultiplicityRatio(inclHist)
        inclHist.Write()
        exclRatio.Write()
        inclRatio.Write()



fileIn = ROOT.TFile.Open("../Unfolding/Zmumu/UnfoldingOutputLep2EtaUnfoldModel.root")
fileOutCtrl = ROOT.TFile.Open("smoothingControlPlots_Lep2EtaUnfoldModel.root", "RECREATE")
fileOut = ROOT.TFile.Open("../Unfolding/Zmumu/UnfoldingOutputLep2EtaUnfoldModel_Smoothed.root", "RECREATE")


sysNames = ["UNFOLDMODEL"]


SetPlotStyle()
copyHistos(fileIn,fileOut)


sys = []

for sysName in sysNames:
    sys.append(GausKernelSmoother.Systematic())
    sys[-1].inputFile = fileIn
    sys[-1].outputControlFile = fileOutCtrl
    sys[-1].outputFile = fileOut
    sys[-1].histoName = "hUnfold"
    sys[-1].variableName = "Lep2Eta"
    sys[-1].sysName = sysName
    sys[-1].logScale = False
    sys[-1].gausWidth = 0.3

tex = LatexDocument.LatexDocument()
tex.texFileName = "smoothedSys_UnfoldModel_Lep2Eta.tex"
tex.figDir = "plots/"
tex.openDocument()

for s in sys:
    print "Smoothing "+s.sysName+"/"+s.variableName
    s.retrieveHistos()
    s.smooth()
    #s.saveControlHistos(tex.figDir)
    s.saveControlHistos(0)
    s.saveNewSys()
    tex.insertFigure(s.sysName+"_"+s.variableName+"_ctrl.eps", s.sysName.replace("_", "-")+", "+s.histoName.replace("_", "-"))

tex.closeDocument()
tex.compile()


fileOutCtrl.Close()
fileOut.Close()
fileIn.Close()
