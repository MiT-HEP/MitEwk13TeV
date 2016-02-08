#!/usr/bin/env python

#import os,sys
from subprocess import call, check_output
import re
import threading
import time
import os,sys

from optparse import OptionParser, OptionGroup

usage = "usage: %prog [options] zmm_eos.conf"
parser=OptionParser(usage=usage)
parser.add_option("-n","--nserver" ,dest='nserver',type='int',help="Number of server to use [Default=%default]",default=10)
parser.add_option("-m","--macro"  ,dest='macro',type='string',help="Macro [Default=%default]",default='selectZee.C')
parser.add_option("-o","--outdir" ,dest='outdir',type='string',help="outdir [Default=%default]",default='Zee')
(opts,args)=parser.parse_args()

def scout():
	''' scout lxplus servers '''
	servers={}
	host = "host lxplus.cern.ch | grep 'has address' | cut -d' ' -f 4 | while read ip ; do host $ip | sed 's/.*pointer //' | sed 's/\.$//'; done "
	hosts = check_output(host, shell=True)
	for l in hosts.split():
		servers[ l ] = 0
	return servers

def PrintServers(servers):
	print "I found the following lxplus:"
	for s in servers:
		print "\t*",s
	print "You want to run on", opts.nserver, "servers and I found ", len(servers)

def MergeDict(dict1,dict2):
	for key in dict2:
		if key not in dict1: dict1[key] = dict2[key]
def TryConnections(servers):
	for s in servers:
		exe="ssh %s mkdir /tmp/%s"%(s,os.environ['USER']) ## the tmp dir may not exists w/o a login
		print "trying: ",exe
		call(exe,shell=True)


servers=scout()

PrintServers(servers)
raw_input("Is ok?")

TryConnections(servers)

def call_exe(exe,server):
	cmd = "ssh " + server + ' "' + exe + '"'
	print " -> Calling: "+ exe
	call(exe,shell=True)
	print " -> Done "+exe
	servers[server] -= 1  ## not completely thread safe if multicore

threads=[]


conf=open( args[0] , "r") 

##cfg = [ ('$ data 1 ... ' , [ file1, file2 ] ) ,
#         ('$ wz ...', [file1,...] )
#       ]
cfg = []

for line in conf:
	line = line.replace('\n','')
	if line[0] == '#' : continue
	if line[0] == '$':
		cfg.append((line,[]))
		continue
	cfg[-1][1].append(line)

def PrintCfg(outname, sampleIdx=0):
	out = open(outname,"w")
	if sampleIdx != 0 : 
		out.write(cfg[0][0] + "\n" ) 

	out.write(cfg[sampleIdx][0] + "\n" ) 

	out.write( '\n'.join( cfg[sampleIdx][1] ) )
	out.close()		
	time.sleep(1) ## make sure file is close, TODO
	return

## while True:
##     try:
##         with open(filename, 'rb') as _:
##             break
##     except IOError:
##         time.sleep(3)

def CompileMacros(macro):
	if True: ##VERBOSE
		print "--> Compiling ", macro
	cmd = "cd " + os.environ['PWD'] + " ; "
	cmd += "eval `scramv1 runtime -sh` ; "  #cmsenv
	cmd += 'root -l <<EOF \n.L %s++\n.q\nEOF'%(macro)
	status = call(cmd,shell=True)
	if status :
		raise Exception("COMPILATION FAILED")
	return
		

def PrepareExe(sampleIdx,server):
	outname= "/tmp/" + os.environ['USER'] + "/%s.%d"%(args[0],sampleIdx)
	PrintCfg(outname,sampleIdx)
	if True:##VERBOSE
	    print " TRANSFERRING CFG '%s' to %s:"%(outname,server)
	    call(["cat",outname])
	    print "\n----------------------------"
	## transfer cfg
	cmd = "rsync -avP %s %s:%s"%(outname,server,outname)
	status = call(cmd,shell=True)
	if status :
		print "Transfer cmd was: '%s'"%cmd
		raise Exception("TRANSFER FAILED")
	##prepare exe	
	exe = "cd " + os.environ['PWD'] + " ; "
	exe += "eval `scramv1 runtime -sh` ; "  #cmsenv
	exe += 'root -l -q \'%s+("%s","%s",0)\' '%(opts.macro,outname,opts.outdir)
	return exe

CompileMacros(opts.macro)
#for exe in args:
for sampleIdx in range(0,len(cfg)):
	server = ""
	## don't send too many job to a single lxplus
	while server == "":
		while threading.activeCount() >= opts.nserver: 
			print "sleep ...",
			time.sleep(3)
			print "wake up"
		for s in servers: 
			if  servers[s] < 1 : 
				server = s
				servers[s] = 1
				break
	#prepare exe line	
	exe= PrepareExe(sampleIdx,server)
	t = threading.Thread(target=call_exe, args= (exe,server)  )
	t.start()
	threads.append( t )

## make sure everything is done
for t in threads:
	t.join()

print " ************* Done ***************"

