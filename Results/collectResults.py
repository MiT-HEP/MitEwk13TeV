#!/bin/env python
import sys,os,re
import math
from subprocess import call,check_output

####### PARSE ARGS ########
from optparse import OptionParser

parser=OptionParser()
parser.add_option("-d","--directory",type="string",help="scouting directory [%default]",default="/afs/cern.ch/work/s/sabrandt/public/SM/Differential/CMSSW_7_6_3_patch2/src/MitEwk13TeV/SignalExtraction")
parser.add_option("-c","--central",type="string",help="central value [%default]",default="Central_Charge")
#parser.add_option("","--s1",type="string",help="systematics one band [%default]",default="Ewk_Free,QCD_Free,Recoil_RooKeys,Recoil_Inclusive")
parser.add_option("","--s1",type="string",help="systematics one band [%default]",default="QCD_?ree,Ewk_Free,Recoil_RooKeys,Recoil_Inclusive")
parser.add_option("","--s2",type="string",help="systematics two band [%default]",default="Pileup")
#parser.add_option("-e","--efficiency",type="string",help="scouting directory [%default]",default="/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Results/%s/ChargeDependentEff/binned.txt")

opts,args=parser.parse_args()

###########################

#cmd="cat /afs/cern.ch/work/s/sabrandt/public/SM/Differential/CMSSW_7_6_3_patch2/src/MitEwk13TeV/SignalExtraction/Wmunu_Central_noCharge/fitresWmm*txt| grep Signal | head -n 1"
def AddSpace(s,n=20):
	toprint=s[:]
	l=len(s)
	for i in range(l,n):
		toprint+=" "
	return toprint


def ProduceTot(d):
	d["Wenu"]["Wtot"] = d["Wenu"]["Wp"]+d["Wenu"]["Wm"]
	d["Wmunu"]["Wtot"] = d["Wmunu"]["Wp"]+d["Wmunu"]["Wm"]
	# ratio
	d["Wenu"]["Wratio"] = d["Wenu"]["Wp"]/d["Wenu"]["Wm"]
	d["Wmunu"]["Wratio"] = d["Wmunu"]["Wp"]/d["Wmunu"]["Wm"]

def ProduceRel(a,b):
	r={}
	for ch in ["Wmunu","Wenu"]:
	   r[ch]={}
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   r[ch][w]=abs(a[ch][w]-b[ch][w])/b[ch][w]
	return r

def ProduceDiff(a,b,c):
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


Central=ReadDict(opts.central)
CentralError=ReadDict(opts.central,True)
#CentralError={}
#for ch in Central:
#	CentralError[ch]={}
#	for w in Central[ch]: 
#		if w=="Wratio":
#			pass
#		else:
#			CentralError[ch][w] = math.sqrt(Central[ch][w])
#	if "Wratio" in Central[ch]:
#		w="Wratio"
#		CentralError[ch][w] = Central[ch][w]* math.sqrt(1./ Central[ch]["Wp"] + 1./ Central[ch]["Wm"] )
for ch in Central:
	if "Wtot" in Central[ch]:
		w="Wtot"
		CentralError[ch][w] = math.sqrt(CentralError[ch]["Wp"]**2 + CentralError[ch]["Wm"]**2 )
	if "Wratio" in Central[ch]:
		w="Wratio"
		CentralError[ch][w] = Central[ch][w]* math.sqrt( (CentralError[ch]["Wp"]/ Central[ch]["Wp"])**2 + (CentralError[ch]["Wm"]/ Central[ch]["Wm"])**2 )

Efficiency=ReadEffXAcc("ChargeDependentEff")
EffErr=ReadEffXAcc("ChargeDependentEff",True)
#
EffPuUp=ReadEffXAcc("PileupUp")
EffPuDown=ReadEffXAcc("PileupDown")
EffPileup=ProduceDiff(EffPuUp,EffPuDown,Efficiency)

Systematics={}
for s in opts.s1.split(","):
	tmp=ReadDict(s)
	Systematics[s] = ProduceRel(tmp,Central)

for s in opts.s2.split(","):
	tmpU=ReadDict(s+"_Up")
	tmpD=ReadDict(s+"_Down")
	Systematics[s] = ProduceDiff(tmpU,tmpD,Central)

def PrintLine(info,d,form="",extra=""):
	print info,
	for ch in ["Wmunu","Wenu"]:
	   for w in ["Wp","Wm","Wtot","Wratio"]:
		   if form=="":
			if d[ch][w]<-99.: print "   ---------",
			else :  print "  ",d[ch][w],
		   else:
			if extra !="" and extra != "x" and w=="Wratio": 
				print "  ",extra%d[ch][w],
			elif extra !="" and extra == "x" and w=="Wratio": 
				print "   ---------",
			elif d[ch][w]<-99.:
				print "   ---------",
			else: print "  ",form%d[ch][w],
	print  ##EOL
	return	



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
for s in Systematics:
	toprint=AddSpace(s,20)
	PrintLine(toprint,Systematics[s],"%10.4f")
	
