#!/usr/bin/env python
import os
import subprocess
import commands

#efftypes = ['Sta','SIT']
#efftypes = ['SIT']
efftypes = ['GsfSel']
for efftype in efftypes:
#    for binnum in range(24):
    for binnum in range(36):
#    for binnum in [0]:
            if binnum == 2:
                continue
            elif binnum == 9:
                continue
            elif binnum == 14:
                continue
            elif binnum == 21:
                continue
            elif binnum == 26:
                continue
            elif binnum == 33:
                continue

            cmd = "/afs/cern.ch/user/x/xniu/WZXSection/CMSSW_7_6_3_patch2/src/MitEwk13TeV/XControl/runShapeUncStep23_bkg.sh Ele "+efftype+" etapt "+str(binnum)+" 1000"
            print cmd
            bsubs_cmd = "bsub -q 8nh -R 'pool > 4000' -C 0 -o" + \
                "/afs/cern.ch/work/x/xniu/public/WZXSection/wz-efficiency/Ele"+efftype+"Eff/ploutputetapt"+str(binnum)+" "+cmd
            status,output=commands.getstatusoutput(bsubs_cmd)
            print output
