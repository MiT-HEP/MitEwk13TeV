"""@package  GausKernelSmoother.py

Smoother: smooth histogram using a Gaussian kernel (based on code from Dag Gillberg)
    * self.histo: histo to be smoothed
    * self.weights: weight for each bin to be used by the Gaussian kernel
    * self.smoothHisto: smoothed histo
    * self.logScale: enable/disable log scale in the Gaussian kernel (log scale is for pT-like variables)
    * self.gausWidth: width of the Gaussian kernel

Systematic: store and compute information for a given systematic source (nominal histo, shifted histo, systematic shift, smoothed systematic shift, shifted histo corresponding to the smoothed systematic shift).
It assumes that the nominal histo and the original shifted histo are stored in a single ROOT file, the nominal histo being in the root directory, and the shifted histo being in a sub-directory (with the same name as the nominal histo)
    * self.inputFile: name of the file containing the nominal histo and the original shifted histo
    * self.outputFile: name of the output file
    * self.outputControlFile: output file containing control plots
    * self.histoName: name of the histo 
    * self.sysName: name of the systematic source (=name of the sub-directory where the corresponding histo is stored)
    * self.histoNom: nominal histo
    * self.histoSys: shifted histo
    * self.histoShift: relative shift between self.histoSys and self.histoNom
    * self.histoWeights: weights to be used for the smoothing (the stat uncert of self.histoNom are used to derive these weights, see createWeightsFromErrors)
    * self.histoSmooth: smoothed relative shift
    * self.histoSmoothContinuous: a fine binned version of self.histoSmooth
    * self.histoSysSmooth: shifted histo computed from self.histoNom and self.histoSmooth

"""



import ROOT
import math
import array

class Smoother:
    def __init__(self):
        self.histo = None
        self.weights = None
        self.smoothHisto = None
        self.logScale = False
        self.gausWidth = 1.

    def getSmoothedValue(self, x):
        if self.logScale and x<=0.:
            raise StandardError("ERROR: use log scale and x<=0.")
        sumw = 0.
        sumwy = 0.
        nbins = self.histo.GetNbinsX()
        for b in range(1,nbins+1):
            xi = self.histo.GetXaxis().GetBinCenter(b)
            yi = self.histo.GetBinContent(b)
            if self.logScale and xi<=0.:
                raise StandardError("ERROR: use log scale and xi<=0.")
            dx = 0.
            if self.logScale:
                dx = (math.log(x) - math.log(xi))/self.gausWidth
            else:
                dx = (x-xi)/self.gausWidth
            wi = ROOT.TMath.Gaus(dx)
            if self.weights:
                wi *= self.weights.GetBinContent(b)
            sumw += wi
            sumwy += wi*yi
        value = 0.
        if sumw>0.:
            value = sumwy/sumw
        return value

    def computeSmoothHisto(self):
        if not self.histo:
            raise StandardError("ERROR: non existing input histo")
        self.smoothHisto = self.histo.Clone(self.histo.GetName()+"_smooth")
        self.smoothHisto.__class__ = ROOT.TH1D
        self.smoothHisto.SetDirectory(0)
        nbins = self.smoothHisto.GetNbinsX()
        for b in range(1,nbins+1):
            x = self.smoothHisto.GetBinCenter(b)
            smoothedValue = self.getSmoothedValue(x)
            self.smoothHisto.SetBinContent(b,smoothedValue)

    def getContinuousSmoothHisto(self):
        if not self.histo:
            raise StandardError("ERROR: non existing input histo")
        mini = self.histo.GetXaxis().GetBinLowEdge(1)
        maxi = self.histo.GetXaxis().GetBinUpEdge(self.histo.GetNbinsX())
        if self.logScale and mini<0.:
            raise StandardError("ERROR: use log scale and min value<0")
        if mini==0.:
            mini = self.histo.GetXaxis().GetBinUpEdge(1)/10.
        nbins = 1000
        bins = []
        if self.logScale:
            dx = (math.log(maxi) - math.log(mini))/nbins
            bins = [math.exp(math.log(mini)+i*dx) for i in range(0,nbins+1)]
        else:
            dx = (maxi - mini)/nbins
            bins = [mini+i*dx for i in range(0,nbins+1)]
        smoothHisto = ROOT.TH1D(self.histo.GetName()+"_cont",self.histo.GetTitle(), nbins, array.array('f',bins))
        for b in range(1,nbins+1):
            x = smoothHisto.GetBinCenter(b)
            smoothedValue = self.getSmoothedValue(x)
            smoothHisto.SetBinContent(b,smoothedValue)
        return smoothHisto


class Systematic:
    def __init__(self):
        self.histoName = ""
        self.sysName = ""
        self.inputFile = None
        self.outputFile = None
        self.outputControlFile = None
        self.histoNom = None
        self.histoSys = None
        self.histoShift = None
        self.histoWeights = None
        self.histoSmooth = None
        self.histoSmoothContinuous = None
        self.histoSysSmooth = None
        self.logScale = False
        self.gausWidth = 0.4

    def retrieveHistos(self):
        if not self.inputFile:
            raise StandardError("ERROR: non existing input file")
        self.histoNom = self.inputFile.Get(self.histoName)
        self.histoNom.__class__ = ROOT.TH1D
        self.histoNom.SetDirectory(0)
        self.histoSys = self.inputFile.Get(self.sysName+"/"+self.histoName)
        self.histoSys.__class__ = ROOT.TH1D
        self.histoSys.SetDirectory(0)
        self.createShiftHisto()
        self.createWeightsFromErrors()

    def createShiftHisto(self):
        if not self.histoName or not self.histoSys:
            raise StandardError("ERROR: non existing histo")
        self.histoShift = self.histoSys.Clone(self.histoSys.GetName()+"_shift")
        self.histoShift.__class__ = ROOT.TH1D
        self.histoShift.SetDirectory(0)
        self.histoShift.Add(self.histoNom,-1.)
        self.histoShift.Divide(self.histoNom)
        for b in range(1,self.histoShift.GetNbinsX()+1):
            if self.histoShift.GetBinContent(b)<0.:
                self.histoShift.SetBinContent(b,self.histoShift.GetBinContent(b)*(-1.))
        for b in range(1,self.histoShift.GetNbinsX()+1):
            self.histoShift.SetBinError(b,0.)

    def createWeightsFromErrors(self):
        if not self.histoName:
            raise StandardError("ERROR: non existing histo")
        self.histoWeights = self.histoNom.Clone(self.sysName+"_"+self.histoNom.GetName()+"_weights")
        self.histoWeights.__class__ = ROOT.TH1D
        self.histoWeights.SetDirectory(0)
        weights = []
        for b in range(1,self.histoNom.GetNbinsX()+1):
            weight = 0.
            if self.histoNom.GetBinContent(b)!=0.:
                relErr = self.histoNom.GetBinError(b)/self.histoNom.GetBinContent(b)
                weight = 1./(relErr**2)
            weights.append(weight)
        sumWeights = sum(weights)
        for b in range(1,self.histoWeights.GetNbinsX()+1):
            self.histoWeights.SetBinContent(b,weights[b-1]/sumWeights)
            self.histoWeights.SetBinError(b,0.)

    def smooth(self):
        smoother = Smoother()
        smoother.logScale = self.logScale
        smoother.gausWidth = self.gausWidth
        smoother.histo = self.histoShift
        smoother.weights = self.histoWeights
        smoother.computeSmoothHisto()
        self.histoSmooth = smoother.smoothHisto
        self.histoSmoothContinuous = smoother.getContinuousSmoothHisto()
        self.createNewSys()

    def createNewSys(self):
        self.histoSysSmooth = self.histoSys.Clone(self.sysName+"_"+self.histoSys.GetName()+"_smooth")
        self.histoSysSmooth.__class__ = ROOT.TH1D
        self.histoSysSmooth.SetDirectory(0)
        nbins = self.histoSysSmooth.GetNbinsX()
        for b in range(1,nbins+1):
            nom = self.histoNom.GetBinContent(b)
            relShiftSmooth = self.histoSmooth.GetBinContent(b)
            sysSmooth = nom*(1.+relShiftSmooth)
            self.histoSysSmooth.SetBinContent(b,sysSmooth)

    def saveControlHistos(self, y):
        if not self.outputControlFile:
            raise StandardError("ERROR: non existing output control file")
        self.outputControlFile.cd()
        canvas = ROOT.TCanvas(self.sysName+"_"+self.variableName+"_ctrl", self.sysName+" shift "+self.histoNom.GetTitle(), 900, 900)
        canvas.SetLogx(y)
        self.histoShift.SetLineWidth(2)
        self.histoShift.SetLineColor(1)
        self.histoShift.GetYaxis().SetTitle("Relative Systematic Shift")
        #self.histoShift.GetXaxis().SetTitle(self.variableName)
        self.histoShift.Draw()
        self.histoSmoothContinuous.SetLineColor(ROOT.kBlue)
        self.histoSmoothContinuous.SetLineWidth(2)
        self.histoSmoothContinuous.SetLineStyle(2)
        self.histoSmoothContinuous.Draw("same")
        self.histoSmooth.SetLineColor(ROOT.kRed)
        self.histoSmooth.SetLineWidth(2)
        self.histoSmooth.Draw("same")
        self.legend=ROOT.TLegend(0.2,0.75,0.35,0.88)
        self.legend.SetTextSize(0.04);
        self.legend.SetTextFont(42);
        self.legend.AddEntry(self.histoShift,"Before Smoothing", "l")
        self.legend.AddEntry(self.histoSmoothContinuous,"Smoothing Result", "l")
        self.legend.AddEntry(self.histoSmooth,"After Smoothing", "l")
        self.legend.Draw()
        canvas.RedrawAxis()
        canvas.Write()
        canvas.Print("plots/"+canvas.GetName()+".png")
        canvas.Print("plots/"+canvas.GetName()+".eps")
        self.histoWeights.Write()

    def saveNewSys(self):
        if not self.outputFile:
            raise StandardError("ERROR: non existing output file")
        self.outputFile.cd()
        if not self.outputFile.Get("SMOOTH_"+self.sysName):
            self.outputFile.mkdir("SMOOTH_"+self.sysName)
        self.outputFile.cd("SMOOTH_"+self.sysName)
        self.histoSysSmooth.SetName(self.histoName)
        self.histoSysSmooth.Write()


