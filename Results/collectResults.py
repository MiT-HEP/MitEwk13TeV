#!/bin/env python
import sys,os,re
import math
from subprocess import call,check_output

####### PARSE ARGS ########
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-d","--directory",type="string",help="scouting directory [%default]",default="/afs/cern.ch/work/s/sabrandt/public/SM/Differential/CMSSW_7_6_3_patch2/src/MitEwk13TeV/SignalExtraction")
parser.add_option("-c","--central",type="string",help="central value [%default]",default="Central_Charge")

opts,args=parser.parse_args()

###########################

def AddSpace(s,n=20):
	toprint=s[:]
	l=len(s)
	for i in range(l,n):
		toprint+=" "
	return toprint

def ProduceTot(d):
	'''Handle and complete dictionaries'''
	d["Wenu"]["Wtot"] = d["Wenu"]["Wp"]+d["Wenu"]["Wm"]
	d["Wmunu"]["Wtot"] = d["Wmunu"]["Wp"]+d["Wmunu"]["Wm"]
	# ratio
	d["Wenu"]["Wratio"] = d["Wenu"]["Wp"]/d["Wenu"]["Wm"]
	d["Wmunu"]["Wratio"] = d["Wmunu"]["Wp"]/d["Wmunu"]["Wm"]

def ProduceRel(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   r[ch][w]=abs(a[ch][w]-b[ch][w])/b[ch][w]
	return r

def ProduceDiffRel(a,b,c):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=max(abs(a[ch][w]-c[ch][w])/(c[ch][w]), abs(b[ch][w]-c[ch][w])/(c[ch][w]))
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceDiff(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]-b[ch][w])
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceProd(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]*b[ch][w])
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceRatio(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]/b[ch][w])
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

def SqrtSum(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=math.sqrt(a[ch][w]**2+b[ch][w]**2)
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

def SqrtSumL(l):
	'''Handle and complete dictionaries'''
	r={}
	for a in l:
	  for ch in ["Wmunu","Wenu"]:
	   if ch not in r : r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   if w not in r[ch] : r[ch][w]=0.;
		   try:
			if a[ch][w]<-99. or r[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=math.sqrt(a[ch][w]**2+r[ch][w]**2)
		   except KeyError: 
			   print "Ch",ch,"w=",w,"not in dictionary"
	return r

## get central value

# Central is here as reference structure, it will be overwritten
Central={"Wmunu":{"Wp":-1,"Wm":-1,"Wtot":-1},"Wenu":{"Wp":-1,"Wm":-1,"Wtot":-1}}
Info= {"Wmunu":{"Wp":"W+(mu)","Wm":"W-(mu)","Wtot":"W(mu)","Wratio":"W+/W-(mu)"},"Wenu":{"Wp":"W+(e)","Wm":"W-(e)","Wtot":"W(e)","Wratio":"W+/W-(e)"}}

for ch in Info:
	for w in Info[ch]:
		Info[ch][w]=AddSpace(Info[ch][w],10)

def ReadDict(select,error=False):
	'''Stephanie'''
	d={}
	for ch in ["Wmunu","Wenu"]:
	   d[ch]={}
	   for w in ["Wp","Wm"]:
	   #for w in ["Wp","Wm","Wtot"]: ## from FIT part1
		fname="fitresW";
		# add channel info
		if ch=="Wmunu": fname += "m"
		else: fname +="e"
	
		# add charge info
		if w == "Wp": fname +="p"
		elif w=="Wm": fname +="m"
		else: pass
		
		fname+='*.txt'
	
		cmd="cat "+opts.directory +"/" + ch + "_" + select + "/"+fname + "| grep Signal | head -n 1 | tr -s ' ' | sed 's/.*://'"
		valString=check_output(cmd,shell=True)
		
		print>>sys.stderr, "DEBUG: cmd='"+cmd+"'"
		print>>sys.stderr, "DEBUG: valString='"+valString+"'"

		if error:
			d[ch][w] = float (valString.split()[2])
		else:
			d[ch][w] = float (valString.split()[0])

	ProduceTot(d) ## from FIT part2/comment
	return d

def ReadEffXAcc(directory="ChargeDependentEff",error=False):
	'''Xinmei'''
	d={}
	for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
		if 'Z' in ch: l = [""]
		if "W" in ch: l = ["Wp","Wm","Wtot"]
		d[ch]={}
		for w in l:
			 #We   Wm   Zee
			 #Wem  Wmm  Zmm
			 #Wep  Wmp
			dirName=""
			if 'Z' in ch : 
				if 'mu' in ch : dirName = "Zmm"
				else: dirName="Zee"
				cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"/"+directory+"/binned.txt | grep 'SF corrected' | sed 's/^.*://'"
			if 'W' in ch : 
				if 'mu' in ch : dirName = "Wm"
				else: dirName="We"

				if w =="Wp": dirName+="p"
				elif w=="Wm": dirName +="m"
				else : pass
				cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"/"+directory+"/sel.txt | grep 'eff corr' | grep total| sed 's/^.*==>//'"

			if dirName == "Wm" or dirName== "We":
				print "FIXME",ch,w
				d[ch][w] = -999.
			else:
				valString=check_output(cmd,shell=True)
				print>>sys.stderr, "DEBUG: cmd='"+cmd+"'"
				print>>sys.stderr, "DEBUG: valString='"+valString+"'"
				if error:
					d[ch][w] = float (valString.split()[2])
				else:
					d[ch][w] = float (valString.split()[0])
	return d

### Read Central value for W
Yields={}
Systematics={}

Central=ReadDict(opts.central)
CentralError=ReadDict(opts.central,True)
Yields["Central"]= Central

## Compute Wtot and errors
for ch in Central:
	if "Wtot" in Central[ch]:
		w="Wtot"
		CentralError[ch][w] = math.sqrt(CentralError[ch]["Wp"]**2 + CentralError[ch]["Wm"]**2 )
	if "Wratio" in Central[ch]:
		w="Wratio"
		CentralError[ch][w] = Central[ch][w]* math.sqrt( (CentralError[ch]["Wp"]/ Central[ch]["Wp"])**2 + (CentralError[ch]["Wm"]/ Central[ch]["Wm"])**2 )

## read Efficiencies 
SystematicsEff={}

Efficiency=ReadEffXAcc("ChargeDependentEff")
EffErr=ReadEffXAcc("ChargeDependentEff",True)
#systematics are done wrt to this value
EfficiencyNoCh=ReadEffXAcc("Final")

#read efficiency systematics
EffPuUp=ReadEffXAcc("PileupUp")
EffPuDown=ReadEffXAcc("PileupDown")
EffPileup=ProduceDiffRel(EffPuUp,EffPuDown,EfficiencyNoCh)
SystematicsEff["Pileup"]=EffPileup

#### read systematics for W

syst1=["QCD_?ree","Ewk_Fix","Recoil_RooKeys","Recoil_Inclusive"]
for s in syst1:
	tmp=ReadDict(s)
	Yields[s]=tmp
	Systematics[s] = ProduceRel(tmp,Central)

syst2=["Pileup"]
for s in syst2:
	tmpU=ReadDict(s+"_Up")
	tmpD=ReadDict(s+"_Down")
	if s=="Pileup":
		Systematics[s] = ProduceDiffRel(tmpU,tmpD,Yields["Recoil_Inclusive"])
		TmpUp=ProduceRatio(tmpU,EffPuUp)
		TmpDown=ProduceRatio(tmpD,EffPuDown)
		TmpYield=ProduceRatio(Yields["Recoil_Inclusive"],Efficiency)
		TotPuCorr = ProduceDiffRel(TmpUp,TmpDown,TmpYield)
	else:
		Systematics[s] = ProduceDiffRel(tmpU,tmpD,Central)

######################################## Z ##################################################

# collect zmm and zee yields
#/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt
#/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt
YieldsZ={}

YieldsZ["ee_reco"]       = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Zee event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["ee_reco_error"] = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Zee event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["ee_ewk"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["ee_ewk_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["ee_top"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["ee_top_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["ee_mc"]           = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Zee expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["ee_mc_error"]     = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zeeyield.txt | sed 's:\.$::'| grep 'Zee expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_reco"]       = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Zmm event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_reco_error"] = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Zmm event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_ewk"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_ewk_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_top"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_top_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_mc"]         = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Zmm expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_mc_error"]   = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmmyield.txt | sed 's:\.$::'| grep 'Zmm expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

#### Z yields and error propagation
YieldsZ["Zee"] = YieldsZ["ee_reco"] - YieldsZ["ee_ewk"] - YieldsZ["ee_top"]
YieldsZ["Zee_error"] = math.sqrt( YieldsZ["ee_reco_error"]**2 + YieldsZ["ee_ewk_error"]**2 + YieldsZ["ee_top_error"]**2)

YieldsZ["Zmumu"] = YieldsZ["mumu_reco"] - YieldsZ["mumu_ewk"] - YieldsZ["mumu_top"]
YieldsZ["Zmumu_error"] = math.sqrt( YieldsZ["mumu_reco_error"]**2 + YieldsZ["mumu_ewk_error"]**2 + YieldsZ["mumu_top_error"]**2)


######################################## PRINT ###############################################
def PrintLine(info,d,form="",extra="",scale=1):
	print info,
	for ch in ["Wmunu","Wenu"]:
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   if form=="":
			if d[ch][w]<-99.: print "   ---------",
			else :  print "  ",d[ch][w]*scale,
		   else:
			if extra !="" and extra != "x" and w=="Wratio": 
				print "  ",extra%(d[ch][w]*scale),
			elif extra !="" and extra == "x" and w=="Wratio": 
				print "   ---------",
			elif w not in d[ch]:
				print "   ---------",
			elif d[ch][w]<-99.:
				print "   ---------",
			else: print "  ",form%(d[ch][w]*scale),
	print  ##EOL
	return	

print '#'
print '# Lumi (/fb)  Error (relative)'
print '#'
print '2.3055         0.027'
print '#'
print "      ",'\t'.join(["Zmumu","Zee"])
print "Yields",'\t'.join(["%.0f"%YieldsZ["Zmumu"],"%.0f"%YieldsZ["Zee"]])
print "Errors",'\t'.join(["%.0f"%YieldsZ["Zmumu_error"],"%.0f"%YieldsZ["Zee_error"]])



print '#'
print '# Lumi (/fb)  Error (relative)'
print '#'
print '2.3055         0.027'
print '#'
PrintLine(AddSpace("#Quantity",20),Info)
print '#'
PrintLine(AddSpace("yield",20),Central,"%10.1f","%10.8f")
PrintLine(AddSpace("yield_err",20),CentralError,"%10.1f","%10.8f")
PrintLine(AddSpace("effxacc",20),Efficiency,"%10.4f","x")
PrintLine(AddSpace("effxacc err",20),EffErr,"%10.4f","x")
print '#'
PrintLine(AddSpace("#Systematics",20),Info)
print '#'
PrintLine(AddSpace("effxacc pu",20),EffPileup,"%10.4f","x")
#PrintLine(AddSpace("effxacc scale",20),EffScale,"%10.4f","x")
#for i in ["Bin","BkgUp","Ori","SigUp"]:
#	PrintLine(AddSpace("effxacc "+i ,20),SystematicsEff[i],"%10.4f","x")

for s in Systematics:
	toprint=AddSpace(s,20)
	PrintLine(toprint,Systematics[s],"%10.4f")

print '# -----------------'
PrintLine(AddSpace("tot pu corr",20),TotPuCorr,"%10.4f","x")

## collapse
print 
print 
print 
print '-----------------------------------'
print '-----------------------------------'
print 


Recoil=SqrtSum(Systematics["Recoil_RooKeys"],Systematics["Recoil_Inclusive"])
Bkg = SqrtSum(Systematics["QCD_?ree"], Systematics["Ewk_Fix"])
Tot=SqrtSumL( [Recoil,Bkg,TotPuCorr,EffErr])

PrintLine(AddSpace("#Quantity",20),Info)
PrintLine(AddSpace("Recoil"),Recoil,"%10.4f","",100.)
PrintLine(AddSpace("Bkg"),Bkg,"%10.4f","",100.)
#PrintLine(AddSpace("Pileup"),Systematics["Pileup"],"%10.4f","",100.)
PrintLine(AddSpace("Pileup"),TotPuCorr,"%10.4f","",100.)
PrintLine(AddSpace("EffErr"),EffErr,"%10.4f","",100.)
print '-----------------------------------'
PrintLine(AddSpace("Tot"),Tot,"%10.4f","",100.)

print '-----------------------------------'
print '-----------------------------------'

print "what: Central Value"
print Central
for what in Yields:
	print '-----------------------------------'
	print "what:",what
	print Yields[what]

	
