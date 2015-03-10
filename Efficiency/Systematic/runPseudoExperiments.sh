       file=$1
  pseudoDir=$2
 pseudoType=$3
  binNumber=$4
 sigPassFit=$5
 bkgPassFit=$6
 sigFailFit=$7
 bkgFailFit=$8
  outputDir=$9
     charge=${10}
 mcfilename=${11}

h=`basename $0`
echo "Script:    $h"
echo "Arguments: $*"

# some basic printing
echo " "; echo "${h}: Show who and where we are";
echo " "
echo " user executing: "`id`;
echo " running on    : "`hostname`;
echo " executing in  : "`pwd`;
echo " submitted from: $HOSTNAME";
echo " ";


# initialize the CMSSW environment
echo " "; echo "${h}: Initialize CMSSW (in $CMSSW_BASE)"; echo " "
workDir=`pwd`
cd   $CMSSW_BASE
eval `scram runtime -sh`

RUNDIR=/home/klawhorn/cms/cmssw/031/CMSSW_5_3_10_patch1/src/EWKAna/Efficiency/Systematic

pwd
cd $workDir
if [ "$workDir" != $RUNDIR ]
then
    cp $RUNDIR/doPseudoFits.C ./
    cp $RUNDIR/rootlogon.C ./
    cp $RUNDIR/*.hh ./
    cp $RUNDIR/*.h ./
    cp $RUNDIR/*.cc ./
fi

id=`id -u`
cp ~/.krb5/x509up_u${id} /tmp/
cp ~/.krb5/krb5cc_${id}  /tmp/krb5cc_${id}
ls -lhrt /tmp/krb5cc_${id}
export KRB5CCNAME=FILE:/tmp/krb5cc_${id}

echo " "; echo "${h}: Starting root job now"; echo " ";

source /scratch/ksung/ROOT/root/bin/thisroot.sh

echo "root -l -q -b rootlogon.C doPseudoFits.C+\(\"$pseudoDir/$pseudoType\",\"$file\",\"$RUNDIR/$pseudoType/analysis/plots/$binNumber.root\",$sigPassFit,$bkgPassFit,$sigFailFit,$bkgFailFit,\"$outputDir\",$charge,\"$mcfilename\"\)"

#root -l -q -b rootlogon.C doPseudoFits.C+\(\"$pseudoDir/$pseudoType\",\"$file\",\"$RUNDIR/$pseudoType/analysis/plots/$binNumber.root\",$sigPassFit,$bkgPassFit,$sigFailFit,$bkgFailFit,\"$outputDir\",$charge,\"$mcfilename\"\)

root -l -q -b rootlogon.C doPseudoFits.C+\(\"$pseudoDir/$pseudoType\",\"$file\",\"$RUNDIR/$pseudoType/analysis/plots/$binNumber.root\",$sigPassFit,$bkgPassFit,$sigFailFit,$bkgFailFit,\"$outputDir\",$charge,\"$mcfilename\"\)

# get the return code from the root job
status=`echo $?`
echo "${h}: Status - $status"

# store the result
#echo " "; echo "${h}: Checking the work area before copy"; echo " "
#ls -lhrt ./
#echo " "; echo "${h}: Checking the remote area before copy (only $dataset file)"; echo " "
#ls -lhrt $outputDir

exit $status


