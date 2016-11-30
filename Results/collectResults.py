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
	try:
		d["Wenu"]["Wratio"] = d["Wenu"]["Wp"]/d["Wenu"]["Wm"]
	except ZeroDivisionError:
		d["Wenu"]["Wratio"] = float(-999)
	try:
		d["Wmunu"]["Wratio"] = d["Wmunu"]["Wp"]/d["Wmunu"]["Wm"]
	except ZeroDivisionError:
		d["Wmunu"]["Wratio"] = float(-999)

def ProduceRel(a,b):
	'''Handle and complete dictionaries''' # ready for Z
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   try:
		   	r[ch][w]=abs(a[ch][w]-b[ch][w])/b[ch][w]
		   except KeyError:
			   pass
	return r

def ProduceDiffRel(a,b,c):
	'''Handle and complete dictionaries''' # ready for Z
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=max(abs(a[ch][w]-c[ch][w])/(c[ch][w]), abs(b[ch][w]-c[ch][w])/(c[ch][w]))
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceDiff(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]-b[ch][w])
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceProd(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]*b[ch][w])
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def ProduceRatio(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=(a[ch][w]/b[ch][w])
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def SqrtSum(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=math.sqrt(a[ch][w]**2+b[ch][w]**2)
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def SqrtSumL(l):
	'''Handle and complete dictionaries'''
	r={}
	for a in l:
	  for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   if ch not in r : r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio",""]:
		   if w not in r[ch] : r[ch][w]=0.;
		   try:
			if a[ch][w]<-99. or r[ch][w]<-99: r[ch][w]=-999.
			else: r[ch][w]=math.sqrt(a[ch][w]**2+r[ch][w]**2)
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

## get central value

# Central is here as reference structure, it will be overwritten
Central={"Wmunu":{"Wp":-1,"Wm":-1,"Wtot":-1},"Wenu":{"Wp":-1,"Wm":-1,"Wtot":-1},"Zmumu":{"":-1},"Zee":{"":-1}}
Info= {"Wmunu":{"Wp":"W+(mu)","Wm":"W-(mu)","Wtot":"W(mu)","Wratio":"W+/W-(mu)"},"Wenu":{"Wp":"W+(e)","Wm":"W-(e)","Wtot":"W(e)","Wratio":"W+/W-(e)"},"Zmumu":{"":"Zmm"},"Zee":{"":"Zee"}}

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
		##except :
		##	print "ERROR-> File doesn't exist ? "
		##	raw_input("ok? I'll try to continue")
		##	valString="0.0 0.0 0.0"
		
		print>>sys.stderr, "DEBUG: cmd='"+cmd+"'"
		print>>sys.stderr, "DEBUG: valString='"+valString+"'"
		
		try:
			if error:
				d[ch][w] = float (valString.split()[2])
			else:
				d[ch][w] = float (valString.split()[0])
		except IndexError:
			print "->ERROR. Try to ignore?"
			raw_input("ok?")
			d[ch][w] = float ( -999)

	ProduceTot(d) ## from FIT part2/comment
	return d

def ReadEffXAcc(directory="ChargeDependentEff",error=False,doOnly=None):
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
			if doOnly != None and (ch,w) not in doOnly : 
				#d[ch][w] = 0.000
				continue

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

				try:
					if error:
						d[ch][w] = float (valString.split()[2])
					else:
						d[ch][w] = float (valString.split()[0])
				except IndexError:
					print "->ERROR. Try to ignore?"
					raw_input("ok?")
					d[ch][w] = float ( -999)
	return d

Yields={}
Systematics={}

######################################## W ##################################################
### Read Central value for W

Central=ReadDict("Central")
CentralError=ReadDict("Central",True)
Yields["Central"]= Central
Yields["CentralError"]= CentralError

## Compute Wtot and errors
for ch in Central:
	if "Wtot" in Central[ch]:
		w="Wtot"
		CentralError[ch][w] = math.sqrt(CentralError[ch]["Wp"]**2 + CentralError[ch]["Wm"]**2 )
	if "Wratio" in Central[ch]:
		w="Wratio"
		CentralError[ch][w] = Central[ch][w]* math.sqrt( (CentralError[ch]["Wp"]/ Central[ch]["Wp"])**2 + (CentralError[ch]["Wm"]/ Central[ch]["Wm"])**2 )

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

YieldsZ["mumu_roch_reco"]       = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Zmm event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_roch_reco_error"] = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Zmm event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_roch_ewk"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_roch_ewk_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'EWK' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_roch_top"]        = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_roch_top_error"]  = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Top' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

YieldsZ["mumu_roch_mc"]         = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Zmm expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[0])
YieldsZ["mumu_roch_mc_error"]   = float(check_output("cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/zmm_yield_rochsys.txt | sed 's:\.$::'| grep 'Zmm expected event yield' | sed 's:^.*is ::' | sed 's:+/-::'",shell=True).split()[1])

#### Z yields and error propagation

Yields["Central"]["Zee"]={}
Yields["Central"]["Zmumu"]={}
Yields["CentralError"]["Zee"]={}
Yields["CentralError"]["Zmumu"]={}

YieldsBkg={}
YieldsBkg["Zee"]={}
YieldsBkg["Zmumu"]={}

Yields["Central"]["Zee"][""] = YieldsZ["ee_reco"] - YieldsZ["ee_ewk"] - YieldsZ["ee_top"]
#Yields["CentralError"]["Zee"][""] = math.sqrt( YieldsZ["ee_reco_error"]**2 + YieldsZ["ee_ewk_error"]**2 + YieldsZ["ee_top_error"]**2)
Yields["CentralError"]["Zee"][""] = YieldsZ["ee_reco_error"] 
YieldsBkg["Zee"][""]= math.sqrt(YieldsZ["ee_ewk_error"]**2 + YieldsZ["ee_top_error"]**2)

Yields["Central"]["Zmumu"][""] = YieldsZ["mumu_reco"] - YieldsZ["mumu_ewk"] - YieldsZ["mumu_top"]
#Yields["CentralError"]["Zmumu"][""] = math.sqrt( YieldsZ["mumu_reco_error"]**2 + YieldsZ["mumu_ewk_error"]**2 + YieldsZ["mumu_top_error"]**2)
Yields["CentralError"]["Zmumu"][""] = YieldsZ["mumu_reco_error"]
YieldsBkg["Zmumu"][""]= math.sqrt(YieldsZ["mumu_ewk_error"]**2 + YieldsZ["mumu_top_error"]**2)

Yields["Roch"]={}
Yields["Roch"]["Zmumu"]={}
Yields["Roch"]["Zmumu"][""] = YieldsZ["mumu_roch_reco"] - YieldsZ["mumu_roch_ewk"] - YieldsZ["mumu_roch_top"]
#Yields["CentralError"]["Zmumu"][""] = math.sqrt( YieldsZ["mumu_reco_error"]**2 + YieldsZ["mumu_ewk_error"]**2 + YieldsZ["mumu_top_error"]**2)
#YieldsBkg["Zmumu"][""]= math.sqrt(YieldsZ["mumu_ewk_error"]**2 + YieldsZ["mumu_top_error"]**2)

#print >>sys.stderr, "XXXX",YieldsBkg["Zmumu"][""],"/", Yields["Central"]["Zmumu"][""]
Systematics["Zbkg_yield"]=ProduceRatio(YieldsBkg, Yields["Central"])
Systematics["Zroch_yield"]=ProduceRel(Yields["Roch"],Yields["Central"])

for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot"]
	for w in l:
		if ch == "Wmunu" or ch== "Wenu":
			Systematics["Zbkg_yield"][ch][w] =0.000

for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Zmumu' in ch: continue
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot"]
	for w in l:
		Systematics["Zroch_yield"][ch][w] =0.000

Systematics["Roch"]={}

for ch in ["Wmunu"]:
	Systematics["Roch"][ch]={}
	for w in ["Wp"]:
		varString=check_output("cat /afs/cern.ch/user/s/sabrandt/work/public/SM/Differential/Results/RochCorrToys/2016_11_30_mu/Pull/pluss.txt | grep mean",shell=True).split()[2]
		Systematics["Roch"][ch][w]=float(varString)
	for w in ["Wm"]:
		varString=check_output("cat /afs/cern.ch/user/s/sabrandt/work/public/SM/Differential/Results/RochCorrToys/2016_11_30_mu/Pull/minus.txt | grep mean",shell=True).split()[2]
		Systematics["Roch"][ch][w]=float(varString)

######################################## SYST EFF ##################################################

## read Efficiencies 
## and compute systematics with respect to that
SystematicsEff={}

Efficiency=ReadEffXAcc("Central")
EffErr=ProduceRatio(ReadEffXAcc("Central",True),Efficiency)
#systematics are done wrt to this value
#EfficiencyNoCh=ReadEffXAcc("Final")

#read efficiency systematics
EffPuUp=ReadEffXAcc("PileupUp")
EffPuDown=ReadEffXAcc("PileupDown")
EffPileup=ProduceDiffRel(EffPuUp,EffPuDown,Efficiency)
SystematicsEff["Pileup"]=EffPileup

EffScaleUp=ReadEffXAcc("ScaleUp")
EffScaleDown=ReadEffXAcc("ScaleDown")
EffScale=ProduceDiffRel(EffScaleUp,EffScaleDown,Efficiency)
SystematicsEff["Scale"]=EffScale

# read bin efficiency
EffBin=ReadEffXAcc("Bin")
SystematicsEff["Bin"]=ProduceRel(EffBin,Efficiency)
###
#EffOri=ReadEffXAcc("Ori")
###
EffSig=ReadEffXAcc("SigUp")
SystematicsEff["Sig"]=ProduceRel(EffSig,Efficiency)
EffBkg=ReadEffXAcc("BkgUp")
SystematicsEff["Bkg"]=ProduceRel(EffBkg,Efficiency)

EffMisId=ReadEffXAcc("CMI",doOnly=[("Wenu","Wp"),("Wenu","Wm")])
#SystematicsEff["CMI"] = ProduceRel(EffMisId,Efficiency)
SystematicsEff["CMI"] = ProduceRel(EffMisId,Efficiency)

for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot"]
	for w in l:
		if (ch,w) not in [("Wenu","Wp"),("Wenu","Wm")]:
			SystematicsEff["CMI"][ch][w] =0.000

#SystematicsEff["CMI"] = ProduceRel(EffMisId,EfficiencyNoCh)

######################################## SYST W ##################################################
#### read systematics for W

syst1=["QCD_?ree","Ewk_Fix","Recoil_RooKeys","Recoil_Inclusive"]
for s in syst1:
	tmp=ReadDict(s)
	Yields[s]=tmp
	Systematics[s] = ProduceRel(tmp,Central)

syst2=["Pileup","Scale"]
for s in syst2:
	if s== "Pileup":
		yieldU=ReadDict(s+"_Up")
		yieldD=ReadDict(s+"_Down")
	else:
		yieldU=ReadDict(s+"up")
		yieldD=ReadDict(s+"down")

	if s=="Pileup":
		for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
			if 'Z' in ch: l = [""]
			if "W" in ch: l = ["Wp","Wm","Wtot"]
			for w in l:
				if ch == "Zmumu" or ch== "Zee":
					yieldU[ch]={}
					yieldU[ch][w]=Yields["Central"][ch][w]
					yieldD[ch]={}
					yieldD[ch][w]=Yields["Central"][ch][w]
					Yields["Recoil_Inclusive"][ch]={}
					Yields["Recoil_Inclusive"][ch][w]=Yields["Central"][ch][w]

		#print "XXX","yieldUp=",yieldU["Zmumu"][""],"effUp",EffPuUp["Zmumu"][""]
		#print "XXX","yieldDn=",yieldD["Zmumu"][""],"effDn",EffPuDown["Zmumu"][""]
		#print "XXX","yield=",Yields["Central"]["Zmumu"][""],"eff",EfficiencyNoCh["Zmumu"][""]

		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Yields["Central"])
		CrossSectionUp=ProduceRatio(yieldU,EffPuUp)
		CrossSectionDown=ProduceRatio(yieldD,EffPuDown)
		CrossSectionYield=ProduceRatio(Yields["Central"],Efficiency)
		#print "XXX","CrossSectionUp=",CrossSectionUp["Zmumu"][""]
		#print "XXX","CrossSectionDn=",CrossSectionDown["Zmumu"][""]
		#print "XXX","CrossSection=",CrossSectionYield["Zmumu"][""]
		TotPuCorr = ProduceDiffRel(CrossSectionUp,CrossSectionDown,CrossSectionYield)

		Wtmp = ProduceDiffRel(yieldU,yieldD,Central)
		for ch in ["Wmunu","Wenu"]:
			if "W" in ch: l = ["Wp","Wm","Wtot"]
			for w in l:
				TotPuCorr[ch][w]=Wtmp[ch][w]
	elif s=="Scale":

		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Yields["Central"])
		CrossSectionUp=ProduceRatio(yieldU,EffScaleUp)
		CrossSectionDown=ProduceRatio(yieldD,EffScaleDown)
		CrossSectionYield=ProduceRatio(Yields["Central"],Efficiency)
		TotScaleCorr = ProduceDiffRel(CrossSectionUp,CrossSectionDown,CrossSectionYield)
		TotScaleCorr["Zee"]={}
		Ztmp= ProduceDiffRel(EffScaleUp,EffScaleDown,Efficiency)
		TotScaleCorr["Zee"][""] = Ztmp["Zee"][""]


	else:
		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Central)


for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot"]
	for w in l:
		if ch == "Zmumu" or ch== "Zee":
			Systematics["QCD_?ree"][ch][w] =0.000

######################################## PRINT ###############################################
def PrintLine(info,d,form="",extra="",scale=1):
	print info,
	totList=[ 
		("Wmunu","Wp"), ("Wmunu","Wm"),
		#("Wmunu","Wtot"),
		("Zmumu",""),
		#("Wmunu","Wratio"),
		("Wenu","Wp"), ("Wenu","Wm"),
		#("Wenu","Wtot"),
		("Zee",""),
		#("Wenu","Wratio"),
		]
	#for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	#   wList=["Wp","Wm","Wtot","Wratio"]
	#   if 'Z' in ch: wList=[""]
	#   for w in wList:
	for ch,w in totList:
		   if form=="":
			if d[ch][w]<-99.: print "   ---------",
			else :  print "  ",d[ch][w]*scale,
		   else:
			if extra !="" and extra != "x" and w=="Wratio": 
				print "  ",extra%(d[ch][w]*scale),
			#elif extra !="" and extra == "x" and w=="Wratio": 
			#	print "   ---------",
			elif ch not in d:
				print "   ---------",
			elif w not in d[ch]:
				print "   ---------",
			elif d[ch][w]<-99.:
				print "   ---------",
			else: print "  ",form%(d[ch][w]*scale),
	print  ##EOL
	return	

############ print '#'
############ print '# Lumi (/fb)  Error (relative)'
############ print '#'
############ print '2.3055         0.027'
############ print '#'
############ print "      ",'\t'.join(["Zmumu","Zee"])
############ print "Yields",'\t'.join(["%.1f"%Yields["Central"]["Zmumu"][""],"%.1f"%Yields["Central"]["Zee"][""]])
############ print "Errors",'\t'.join(["%.1f"%Yields["CentralError"]["Zmumu"][""],"%.1f"%Yields["CentralError"]["Zee"][""]])
############ print "EffXAc",'\t'.join(["%.4f"%Efficiency["Zmumu"][""],"%.4f"%Efficiency["Zee"][""]])
############ print "EffErr",'\t'.join(["%.4f"%EffErr["Zmumu"][""],"%.4f"%EffErr["Zee"][""]])



print '#'
print '# Lumi (/fb)  Error (relative)'
print '#'
print '2.3055         0.027'
print '#'
PrintLine(AddSpace("#Quantity",20),Info)
print '#'
PrintLine(AddSpace("yield",20),Central,"%10.1f","%10.8f")
PrintLine(AddSpace("yield_err",20),CentralError,"%10.1f","%10.8f")
PrintLine(AddSpace("effxacc",20),Efficiency,"%10.4f")
PrintLine(AddSpace("effxacc err",20),EffErr,"%10.4f")
print '#'
PrintLine(AddSpace("#Systematics",20),Info)
print '#'
PrintLine(AddSpace("effxacc pu",20),SystematicsEff["Pileup"],"%10.4f")
PrintLine(AddSpace("effxacc bin",20),SystematicsEff["Bin"],"%10.4f")
PrintLine(AddSpace("effxacc Sig",20),SystematicsEff["Sig"],"%10.4f")
PrintLine(AddSpace("effxacc Bkg",20),SystematicsEff["Bkg"],"%10.4f")
PrintLine(AddSpace("effxacc chmisid",20),SystematicsEff["CMI"],"%10.4f")
PrintLine(AddSpace("effxacc scale",20),SystematicsEff["Scale"],"%10.4f")

#PrintLine(AddSpace("effxacc scale",20),EffScale,"%10.4f","x")
#for i in ["Bin","BkgUp","Ori","SigUp"]:
#	PrintLine(AddSpace("effxacc "+i ,20),SystematicsEff[i],"%10.4f","x")

for s in Systematics:
	toprint=AddSpace(s,20)
	PrintLine(toprint,Systematics[s],"%10.4f")


## collapse
print 
print 
print 
print '-----------------------------------'
print '|  going to print table for plots |'
print '-----------------------------------'
print 



PrintLine(AddSpace("#Quantity",20),Info)

#PrintLine(AddSpace("eff_stat",20),EffErr,"%10.4f")
#PrintLine(AddSpace("lepton_eff",20),SqrtSumL([SystematicsEff["Sig"],SystematicsEff["Bkg"]]),"%10.4f")
#PrintLine(AddSpace("binning_eff",20),SystematicsEff["Bin"],"%10.4f")
#PrintLine(AddSpace("charge_mis",20),SystematicsEff["CMI"],"%10.4f")
#PrintLine(AddSpace("bkg_model"),SqrtSum(Systematics["QCD_?ree"], Systematics["Zbkg_yield"]),"%10.4f")
#PrintLine(AddSpace("ewk_norm"),Systematics["Ewk_Fix"],"%10.4f")
PrintLine(AddSpace("pu_model"),TotPuCorr,"%10.4f","")
PrintLine(AddSpace("scale_electron"),TotScaleCorr,"%10.4f","")
#PrintLine(AddSpace("met_model"),Systematics["Recoil_Inclusive"],"%10.4f","")
#PrintLine(AddSpace("recoil_model"),Systematics["Recoil_RooKeys"],"%10.4f","")




print 
print 
print 
print '-----------------------------------'
print '|  going to print table for latex |'
print '|              %                  |'
print '-----------------------------------'
print 
Recoil=SqrtSum(Systematics["Recoil_RooKeys"],Systematics["Recoil_Inclusive"])
Bkg = SqrtSumL([Systematics["QCD_?ree"], Systematics["Ewk_Fix"],Systematics["Zbkg_yield"]])
Eff = SqrtSumL( [EffErr,SystematicsEff["Sig"],SystematicsEff["Bkg"],SystematicsEff["Bin"],SystematicsEff["CMI"] ]) 
Tot = SqrtSumL( [ Recoil,Bkg,Eff,TotPuCorr ,TotScaleCorr] ) 

PrintLine(AddSpace("effxacc"),Eff,"%10.4f","",100)
PrintLine(AddSpace("bkg_model"),Bkg,"%10.4f","",100)
PrintLine(AddSpace("recoil"),Recoil,"%10.4f","",100)
PrintLine(AddSpace("pu_model"),TotPuCorr,"%10.4f","",100)
PrintLine(AddSpace("scale_electron"),TotScaleCorr,"%10.4f","",100)
print '-----------------------------------'
PrintLine(AddSpace("tot"),Tot,"%10.4f","",100)

	
print '-----------------------------------'
print '-----------------------------------'
