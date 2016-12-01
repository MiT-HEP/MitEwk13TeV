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

def AddSpace(s,n=20,begin=False):
	toprint=s[:]
	l=len(s)
	toadd=""
	for i in range(l,n):
	    toadd += " "

	if begin: toprint = toadd + toprint
	else : toprint = toprint + toadd

	return toprint

def Produce(d, what, check, func):
	''' standard set of checks'''
	mych,myw= what
	if mych not in d : d[mych]={}
	r=0.0
	tofill=True
	for ch,w in check:
		if   ch not in d or \
			w not in d[ch] or\
			d[ch][w]<-99:
			tofill=False
			break
	if tofill:
		try:
			r= func(d)
		except KeyError, ZeroDivisionError:
			r=0.0
	d[mych][myw]=r

############ SINGLE OPERATIONS ##################

def ProduceTotAsAverage(d):
	'''Handle and complete dictionaries'''
	Produce(d, ("Wenu","Wtot"), [("Wenu","Wp"),("Wenu","Wm")],lambda d: (d["Wenu"]["Wp"]+d["Wenu"]["Wm"])/2.)
	Produce(d, ("Wmunu","Wtot"), [("Wmunu","Wp"),("Wmunu","Wm")],lambda d: (d["Wmunu"]["Wp"]+d["Wmunu"]["Wm"])/2.)

def ProduceTotAsSum(d):
	Produce(d, ("Wenu","Wtot"), [("Wenu","Wp"),("Wenu","Wm")] ,lambda d: d["Wenu"]["Wp"]+d["Wenu"]["Wm"])
	Produce(d, ("Wmunu","Wtot"), [("Wmunu","Wp"),("Wmunu","Wm")] ,lambda d: d["Wmunu"]["Wp"]+d["Wmunu"]["Wm"])

def ProduceRatioAsRatio(d):
	'''Handle and complete dictionaries'''
	# ratio
	Produce(d, ("Wenu","Wratio"), [("Wenu","Wp"),("Wenu","Wm")] ,lambda d: d["Wenu"]["Wp"]/d["Wenu"]["Wm"])
	Produce(d, ("Wmunu","Wratio"), [("Wmunu","Wp"),("Wmunu","Wm")] ,lambda d: d["Wmunu"]["Wp"]/d["Wmunu"]["Wm"])

	for ch in ["e","mu"]:
	   if "W"+ch+"nu" not in  d :  d["W"+ch+"nu"]={}
	   for w in ["Wp","Wm","W"]:
		w2=w
		if w=="W": w2="Wtot"
		Produce( d, ( "W"+ch+"nu",w+"OverZ"), [("W"+ch+"nu",w2),("Z"+ch+ch,"")],
				lambda d: d["W"+ch+"nu"][w2]/d["Z"+ch+ch][""]
				)

def ProduceRatioAsSqrtSum(d):
	'''Handle and complete dictionaries'''
	# ratio
	Produce(d, ("Wenu","Wratio"), [("Wenu","Wp"),("Wenu","Wm")] ,lambda d: math.sqrt(d["Wenu"]["Wp"]**2 + d["Wenu"]["Wm"]**2))
	Produce(d, ("Wmunu","Wratio"), [("Wmunu","Wp"),("Wmunu","Wm")] ,lambda d: math.sqrt(d["Wmunu"]["Wp"]**2+d["Wmunu"]["Wm"]**2))
	for ch in ["e","mu"]:
	   if "W"+ch+"nu" not in  d :  d["W"+ch+"nu"]={}
	   for w in ["Wp","Wm","W"]:
		w2=w
		if w=="W": w2="Wtot"
		Produce( d, ( "W"+ch+"nu",w+"OverZ"), [("W"+ch+"nu",w2),("Z"+ch+ch,"")],
				lambda d: math.sqrt(d["W"+ch+"nu"][w2]**2 + d["Z"+ch+ch][""]**2)
				)

def ProduceRatioAsSum(d):
	'''Handle and complete dictionaries'''
	# ratio
	Produce(d, ("Wenu","Wratio"), [("Wenu","Wp"),("Wenu","Wm")] ,lambda d: abs(d["Wenu"]["Wp"] + d["Wenu"]["Wm"]))
	Produce(d, ("Wmunu","Wratio"), [("Wmunu","Wp"),("Wmunu","Wm")] ,lambda d: abs(d["Wmunu"]["Wp"] + d["Wmunu"]["Wm"]))

	for ch in ["e","mu"]:
	   if "W"+ch+"nu" not in  d :  d["W"+ch+"nu"]={}
	   for w in ["Wp","Wm","W"]:
		w2=w
		if w=="W": w2="Wtot"
		Produce( d, ( "W"+ch+"nu",w+"OverZ"), [("W"+ch+"nu",w2),("Z"+ch+ch,"")],
				lambda d: d["W"+ch+"nu"][w2] + d["Z"+ch+ch][""]
				)

def ProduceRatioAsDiff(d):
	'''Handle and complete dictionaries'''
	# ratio
	Produce(d, ("Wenu","Wratio"), [("Wenu","Wp"),("Wenu","Wm")] ,lambda d: abs(d["Wenu"]["Wp"] - d["Wenu"]["Wm"]))
	Produce(d, ("Wmunu","Wratio"), [("Wmunu","Wp"),("Wmunu","Wm")] ,lambda d: abs(d["Wmunu"]["Wp"] -d["Wmunu"]["Wm"]))

	for ch in ["e","mu"]:
	   if "W"+ch+"nu" not in  d :  d["W"+ch+"nu"]={}
	   for w in ["Wp","Wm","W"]:
		w2=w
		if w=="W": w2="Wtot"
		Produce( d, ( "W"+ch+"nu",w+"OverZ"), [("W"+ch+"nu",w2),("Z"+ch+ch,"")],
				lambda d: abs(d["W"+ch+"nu"][w2]-d["Z"+ch+ch][""])
				)
def ZeroNotPresent(d):
	for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
		if 'Z' in ch: l = [""]
		if "W" in ch: l = ["Wp","Wm","Wtot","Wratio","WOverZ","WpOverZ","WmOverZ"]
		if ch not in d: d[ch]={}
		for w in l:
			if w not in d[ch]: d[ch][w]=0.0
############ BINARY OPERATIONS #################

def ProduceRel(a,b):
	'''Handle and complete dictionaries''' # ready for Z
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   try:
		   	r[ch][w]=abs(a[ch][w]-b[ch][w])/b[ch][w]
		   except KeyError:
			   pass
		   except ZeroDivisionError:
		   	r[ch][w]=0.000
	return r

def ProduceDiffRel(a,b,c,centralException=False):
	'''Handle and complete dictionaries''' # ready for Z
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=0.00
			else: r[ch][w]=max(abs(a[ch][w]-c[ch][w])/(c[ch][w]), abs(b[ch][w]-c[ch][w])/(c[ch][w]))
		   except KeyError: 
			   pass
		   except ZeroDivisionError:
			if centralException:
		   		r[ch][w]=c[ch][w]
			else:
		   		r[ch][w]=0.000
	return r

def ProduceDiff(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=0.00
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
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=0.00
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
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   #r[ch][w]=abs(a[ch][w]-b[ch][w])/(c[ch][w]*2.)
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=0.00
			else: r[ch][w]=(a[ch][w]/b[ch][w])
		   except KeyError: 
			   pass
		   except ZeroDivisionError:
			r[ch][w]=0.00
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

def SqrtSum(a,b):
	'''Handle and complete dictionaries'''
	r={}
	for ch in ["Wmunu","Wenu","Zee","Zmumu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   try:
			if a[ch][w]<-99. or b[ch][w]<-99: r[ch][w]=0.00
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
	   for w in ["Wp","Wm","Wtot","Wratio","","WpOverZ","WmOverZ","WOverZ"]:
		   if w not in r[ch] : r[ch][w]=0.;
		   try:
			if a[ch][w]<-99. or r[ch][w]<-99: r[ch][w]=0.00
			else: r[ch][w]=math.sqrt(a[ch][w]**2+r[ch][w]**2)
		   except KeyError: 
			   pass
			   #print "Ch",ch,"w=",w,"not in dictionary"
	return r

## get central value

# Central is here as reference structure, it will be overwritten
Central={"Wmunu":{"Wp":-1,"Wm":-1,"Wtot":-1},"Wenu":{"Wp":-1,"Wm":-1,"Wtot":-1},"Zmumu":{"":-1},"Zee":{"":-1}}
Info= {"Wmunu":{"Wp":"W+(mu)","Wm":"W-(mu)","Wtot":"W(mu)","Wratio":"W+/W-(mu)","WpOverZ":"W+/Z(mu)","WmOverZ":"W-/Z(mu)","WOverZ":"W/Z(mu)"},"Wenu":{"Wp":"W+(e)","Wm":"W-(e)","Wtot":"W(e)","Wratio":"W+/W-(e)","WpOverZ":"W+/Z(e)","WmOverZ":"W-/Z(e)","WOverZ":"W/Z(e)"},"Zmumu":{"":"Zmm"},"Zee":{"":"Zee"}}

#for ch in Info:
#	for w in Info[ch]:
#		Info[ch][w]=AddSpace(Info[ch][w],10)

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

	ProduceTotAsSum(d) ## from FIT part2/comment
	ProduceRatioAsRatio(d)

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

			#if dirName == "Wm" or dirName== "We":
			if dirName== "We": ###FIXME
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

def ReadEffXAccForTot(d,directory="ChargeDependentEff",error=False,doOnly=None):
	'''Xinmei
	-> add the tot eff x acc for W, using W+ and W-
	'''
	for ch in ["Wmunu","Wenu"]:
		if "W" in ch: l = ["Wtot"]
		if ch not in d: d[ch]={}
		for w in l:
			if doOnly != None and (ch,w) not in doOnly : 
				#d[ch][w] = 0.000
				continue

			dirName=""

			if 'mu' in ch : dirName = "Wm"
			else: dirName="We"

			cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"p/"+directory+"/sel.txt | grep 'eff corr' | grep total| sed 's/^.*==>//'"
			eff1s=check_output(cmd,shell=True)

			cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"m/"+directory+"/sel.txt | grep 'eff corr' | grep total| sed 's/^.*==>//'"
			eff2s=check_output(cmd,shell=True)

			print>>sys.stderr, "DEBUG: cmd='"+cmd+"'"
			print>>sys.stderr, "DEBUG: valString='"+eff1s+"'", "valString='"+eff2s+"'"

			cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"p/"+directory+"/sel.txt | grep 'eff corr' | grep total| sed 's/^.*:\ *//g'"
			all1s=check_output(cmd,shell=True)
			cmd="cat /afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/"+dirName+"m/"+directory+"/sel.txt | grep 'eff corr' | grep total| sed 's/^.*:\ *//g'"
			all2s=check_output(cmd,shell=True)


			try:
				all1 = float (all1s.split()[2])
				all2 = float (all2s.split()[2])
				if error:
					eff1 = float (eff1s.split()[2])
					eff2 = float (eff2s.split()[2])
				else:
					eff1 = float (eff1s.split()[0])
					eff2 = float (eff2s.split()[0])
				d[ch][w] = (eff1*all1+eff2*all2)/(all1+all2)
			except IndexError:
				print "->ERROR. Try to ignore?"
				raw_input("ok?")
				d[ch][w] = float ( -999)
				return 
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
ZeroNotPresent(Systematics["Zbkg_yield"])

for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Zmumu' in ch: continue
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot","Wratio"]
	for w in l:
		Systematics["Zroch_yield"][ch][w] =0.000
	if "W" in ch: l = ["WOverZ","WmOverZ","WpOverZ"]
	for w in l:
		Systematics["Zroch_yield"][ch][w]=Systematics["Zroch_yield"]["Zmumu"][""]

Systematics["Roch"]={}

for ch in ["Wmunu"]:
	Systematics["Roch"][ch]={}
	for w in ["Wp"]:
		varString=check_output("cat /afs/cern.ch/user/s/sabrandt/work/public/SM/Differential/Results/RochCorrToys/2016_11_30_mu/Pull/pluss.txt | grep mean",shell=True).split()[2]
		Systematics["Roch"][ch][w]=float(varString)
	for w in ["Wm"]:
		varString=check_output("cat /afs/cern.ch/user/s/sabrandt/work/public/SM/Differential/Results/RochCorrToys/2016_11_30_mu/Pull/minus.txt | grep mean",shell=True).split()[2]
		Systematics["Roch"][ch][w]=float(varString)

ZeroNotPresent(Systematics["Roch"])
######################################## SYST EFF ##################################################

## read Efficiencies 
## and compute systematics with respect to that
SystematicsEff={}

Efficiency=ReadEffXAcc("Central")
ReadEffXAccForTot(Efficiency,"Central",error=False,doOnly=[("Wenu","Wtot")]) ## compute effxacc for Wtot enu, starting form W+ and W-
ProduceRatioAsRatio(Efficiency) ## ratio only

BareEffErr=ReadEffXAcc("Central",True)
ReadEffXAccForTot(BareEffErr,"Central",error=True,doOnly=[("Wenu","Wtot")]) ## compute effxacc for Wtot enu, starting form W+ and W-
EffErr=ProduceRatio(BareEffErr,Efficiency)
#systematics are done wrt to this value
ProduceRatioAsSqrtSum(EffErr)

#read efficiency systematics
EffPuUp=ReadEffXAcc("PileupUp")
EffPuDown=ReadEffXAcc("PileupDown")
EffPileup=ProduceDiffRel(EffPuUp,EffPuDown,Efficiency)
ProduceTotAsAverage(EffPileup)
ProduceRatioAsDiff(EffPileup)
SystematicsEff["Pileup"]=EffPileup

EffScaleUp=ReadEffXAcc("ScaleUp")
EffScaleDown=ReadEffXAcc("ScaleDown")
EffScale=ProduceDiffRel(EffScaleUp,EffScaleDown,Efficiency)
ProduceTotAsAverage(EffScale)
ProduceRatioAsDiff(EffScale)
SystematicsEff["Scale"]=EffScale

for ch in ["Wmunu","Zmumu"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot","Wratio","WpOverZ","WmOverZ","WOverZ"]
	for w in l:
		SystematicsEff["Scale"][ch][w] =0.000

# read bin efficiency
EffBin=ReadEffXAcc("Bin")
SystematicsEff["Bin"]=ProduceRel(EffBin,Efficiency)
ProduceTotAsAverage(SystematicsEff["Bin"]) ##FIXME ???
ProduceRatioAsDiff(SystematicsEff["Bin"])
###
#EffOri=ReadEffXAcc("Ori")
###
EffSig=ReadEffXAcc("SigUp")
SystematicsEff["Sig"]=ProduceRel(EffSig,Efficiency)
ProduceTotAsAverage(SystematicsEff["Sig"])
ProduceRatioAsDiff(SystematicsEff["Sig"])
## bkg
EffBkg=ReadEffXAcc("BkgUp")
SystematicsEff["Bkg"]=ProduceRel(EffBkg,Efficiency)
ProduceTotAsAverage(SystematicsEff["Bkg"])
ProduceRatioAsDiff(SystematicsEff["Bkg"])

EffMisId=ReadEffXAcc("CMI",doOnly=[("Wenu","Wp"),("Wenu","Wm")])
#SystematicsEff["CMI"] = ProduceRel(EffMisId,Efficiency)
SystematicsEff["CMI"] = ProduceRel(EffMisId,Efficiency)
ProduceRatioAsSum(SystematicsEff["CMI"])

for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot","Wratio"]
	for w in l:
		if (ch,w) not in [("Wenu","Wp"),("Wenu","Wm"),("Wenu","Wratio")]:
			SystematicsEff["CMI"][ch][w] =0.000

#SystematicsEff["CMI"] = ProduceRel(EffMisId,EfficiencyNoCh)

######################################## SYST W ##################################################
#### read systematics for W

syst1=["QCD_?ree","Ewk_Fix","Recoil_RooKeys","Recoil_Inclusive"]
for s in syst1:
	tmp=ReadDict(s)
	Yields[s]=tmp
	Systematics[s] = ProduceRel(tmp,Central)
	ProduceRatioAsDiff(Systematics[s])
	ZeroNotPresent(Systematics[s])

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


		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Yields["Central"])
		CrossSectionUp=ProduceRatio(yieldU,EffPuUp)
		CrossSectionDown=ProduceRatio(yieldD,EffPuDown)
		CrossSectionYield=ProduceRatio(Yields["Central"],Efficiency)

		TotPuCorr = ProduceDiffRel(CrossSectionUp,CrossSectionDown,CrossSectionYield)

		Wtmp = ProduceDiffRel(yieldU,yieldD,Central)
		for ch in ["Wmunu","Wenu"]:
			if "W" in ch: l = ["Wp","Wm","Wtot"]
			for w in l:
				TotPuCorr[ch][w]=Wtmp[ch][w]
		ZeroNotPresent(Systematics[s])
		ZeroNotPresent(TotPuCorr)
	elif s=="Scale":

		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Yields["Central"])

		for ch in ["Wmunu","Zee","Zmumu"]:
			if "W" in ch: l = ["Wp","Wm","Wtot","Wratio","WpOverZ","WmOverZ","WOverZ"]
			if "Z" in ch: l = [""]
			for w in l:
				Systematics[s][ch][w]=0.
			if "W" in ch: l = ["Wp","Wm","Wtot"]
			for w in l: ## for the ratios
				yieldU = Yields["Central"]
				yieldD = Yields["Central"]

		ProduceTotAsSum(yieldU)
		ProduceTotAsSum(yieldD)
		ProduceRatioAsRatio(yieldU)
		ProduceRatioAsRatio(yieldD)


		CrossSectionUp=ProduceRatio(yieldU,EffScaleUp)
		CrossSectionDown=ProduceRatio(yieldD,EffScaleDown)
		CrossSectionYield=ProduceRatio(Yields["Central"],Efficiency)
		TotScaleCorr = ProduceDiffRel(CrossSectionUp,CrossSectionDown,CrossSectionYield)
		TotScaleCorr["Zee"]={}
		#Ztmp= ProduceDiffRel(EffScaleUp,EffScaleDown,Efficiency)
		TotScaleCorr["Zee"][""] = SystematicsEff["Scale"]["Zee"][""]
		TotScaleCorr["Wenu"]["Wtot"] = max( SystematicsEff["Scale"]["Wenu"]["Wtot"],Systematics["Scale"]["Wenu"]["Wtot"]) ## FIXME
		ProduceRatioAsDiff(TotScaleCorr)
		ZeroNotPresent(TotScaleCorr)
		ZeroNotPresent(Systematics[s])

		## make sure for mu is 0
		for ch in ["Wmunu","Zmumu"]:
			if "W" in ch: l = ["Wp","Wm","Wtot","Wratio","WpOverZ","WmOverZ","WOverZ"]
			if "Z" in ch: l= [""]
			for w in l:
				TotScaleCorr[ch][w]=0.
	else:
		Systematics[s] = ProduceDiffRel(yieldU,yieldD,Central)
		ProduceRatioAsDiff(Systematics[s])
		ZeroNotPresent(Systematics[s])


for ch in ["Wmunu","Wenu","Zmumu","Zee"]:
	if 'Z' in ch: l = [""]
	if "W" in ch: l = ["Wp","Wm","Wtot"]
	for w in l:
		if ch == "Zmumu" or ch== "Zee":
			Systematics["QCD_?ree"][ch][w] =0.000

######################################## PRINT ###############################################
def PrintLine(info,d,scale=1,bigNumbers=False):
	print info,

	if bigNumbers:
		totList=[ 
		("Wmunu","Wp"), ("Wmunu","Wm"),
		("Wmunu","Wtot"),
		("Zmumu",""),
		("Wenu","Wp"), ("Wenu","Wm"),
		("Wenu","Wtot"),
		("Zee",""),
		]
	else:
		totList=[ 
		("Wmunu","Wp"), ("Wmunu","Wm"),
		("Wmunu","Wtot"),
		("Zmumu",""),
		("Wmunu","Wratio"),
		("Wmunu","WOverZ"),
		("Wenu","Wp"), ("Wenu","Wm"),
		("Wenu","Wtot"),
		("Zee",""),
		("Wenu","Wratio"),
		("Wenu","WOverZ"),
		("Wenu","WpOverZ"),
		("Wenu","WmOverZ"),
		("Wmunu","WpOverZ"),
		("Wmunu","WmOverZ"),
		]
	if bigNumbers:
		form="%10.4f"
		n=10
	else: 
		form="%6.4f"
		n=6

	#extra="%10.8f"
	blank=AddSpace("----",n,True)
	for ch,w in totList:
		if ch not in d: print blank,
		elif w not in d[ch]:  print blank,
		elif d[ch][w]<-99.: print blank,
		elif isinstance(d[ch][w], basestring): 
			print AddSpace(d[ch][w],n,True),
		else: print form%(d[ch][w]*scale),
	print  ##EOL
	return	


## PRODUCE RATIOS


print '#'
print '# Lumi (/fb)  Error (relative)'
print '#'
print '2.3055         0.027'
print '#'
PrintLine(AddSpace("#Quantity",20),Info,1.,True)
print '#'
PrintLine(AddSpace("yield",20),Central,1,True)
PrintLine(AddSpace("yield_err",20),CentralError,1.,True)
PrintLine(AddSpace("effxacc",20),Efficiency,1.,True)
PrintLine(AddSpace("effxacc err",20),EffErr,1.,True)
print '#'
PrintLine(AddSpace("#Systematics",20),Info)
print '#'
PrintLine(AddSpace("effxacc pu",20),SystematicsEff["Pileup"])
PrintLine(AddSpace("effxacc bin",20),SystematicsEff["Bin"])
PrintLine(AddSpace("effxacc Sig",20),SystematicsEff["Sig"])
PrintLine(AddSpace("effxacc Bkg",20),SystematicsEff["Bkg"])
PrintLine(AddSpace("effxacc chmisid",20),SystematicsEff["CMI"])
PrintLine(AddSpace("effxacc scale",20),SystematicsEff["Scale"])


for s in Systematics:
	toprint=AddSpace(s,20)
	PrintLine(toprint,Systematics[s])


## collapse
print 
print 
print 
print '-----------------------------------'
print '|  going to print table for plots |'
print '-----------------------------------'
print 



PrintLine(AddSpace("#Quantity",20),Info)

PrintLine(AddSpace("pu_model"),TotPuCorr)
PrintLine(AddSpace("scale_electron"),TotScaleCorr)



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

PrintLine(AddSpace("effxacc"),Eff,100)
PrintLine(AddSpace("bkg_model"),Bkg,100)
PrintLine(AddSpace("recoil"),Recoil,100)
PrintLine(AddSpace("pu_model"),TotPuCorr,100)
PrintLine(AddSpace("scale_electron"),TotScaleCorr,100)
print '-----------------------------------'
PrintLine(AddSpace("tot"),Tot,100)

	
print '-----------------------------------'
print '-----------------------------------'
