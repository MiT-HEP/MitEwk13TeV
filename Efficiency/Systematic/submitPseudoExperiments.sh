echo " "
echo "Starting data processing with arguments:"
echo "  --> $*"
  pseudoDir=$1
 pseudoType=$2
  binNumber=$3
 sigPassFit=$4
 bkgPassFit=$5
 sigFailFit=$6
 bkgFailFit=$7
  outputDir=$8
     charge=$9
 mcfilename=${10}
#
# ./submitPseudoExperiments.sh /scratch/klawhorn/EffSysStore CB_MuSelEff etapt_6 2 1 2 1 /scratch/klawhorn/EffSysStore/CB_MuSelEffResults 0 /scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probes.root
#

jobId=`date +%j%m%d%k%M%S`

workDir=/home/$USER/cms/condor
RUNDIR=/home/klawhorn/cms/cmssw/031/CMSSW_5_3_10_patch1/src/EWKAna/Efficiency/Systematic
mkdir -p $workDir
cd $workDir
cp $RUNDIR/runPseudoExperiments.sh ./

id=`id -u`
mkdir             -p  ~/.krb5/
cp /tmp/x509up_u${id} ~/.krb5/
KRB5CCNAME=`klist -5 | grep 'Ticket cache:' | cut -d' ' -f 3`
if ! [ -z $KRB5CCNAME ]
then
  mkdir    -p  ~/.krb5/
  chmod 0      ~/.krb5
  chmod u=rwx  ~/.krb5
  file=`echo $KRB5CCNAME | cut -d: -f2`
  if [ -f "$file" ]
  then
    cp $file ~/.krb5/krb5cc_${id}
  else
    echo " ERROR -- missing kerberos ticket ($KRB5CCNAME)."
    exit 1
  fi
else
  echo " ERROR -- missing kerberos ticket ($KRB5CCNAME)."
  exit 1
fi

echo
echo

for file in `ls $pseudoDir/$pseudoType/ | grep ^${binNumber}_`
do

  echo $file

  logFile=`echo $file | tr '/' '+'`
  logFile=/tmp/$USER/$logFile
  mkdir -p /tmp/$USER
  rm    -f $logFile

  echo "runPseudoExperiments.sh $file $pseudoDir $pseudoType $binNumber $sigPassFit $bkgPassFit $sigFailFit $bkgFailFit $outputDir $charge $mcfilename  "

  #./runPseudoExperiments.sh $file $pseudoDir $pseudoType $binNumber $sigPassFit $bkgPassFit $sigFailFit $bkgFailFit $outputDir $charge $mcfilename

cat > submit.cmd <<EOF
Universe                = vanilla
Requirements            = ((Arch == "X86_64") && (Machine != "t3btch112.mit.edu") && (Disk >= DiskUsage) && ((Memory * 1024) >= ImageSize) && (HasFileTransfer))
Notification            = Error
Executable              = runPseudoExperiments.sh
Arguments               = $file $pseudoDir $pseudoType $binNumber $sigPassFit $bkgPassFit $sigFailFit $bkgFailFit $outputDir $charge $mcfilename
Rank                    = Mips
GetEnv                  = True
Initialdir              = $workDir
Input                   = /dev/null
Output                  = ${outputDir}/${file}.out
Error                   = ${outputDir}/${file}.err
Log                     = $logFile
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT
Queue
EOF
  #echo "before"
  condor_submit submit.cmd >& /dev/null;
  #echo "after"
  rm submit.cmd

done

exit 0

# actually do pseudofits                                                                                                                                                                                   
#root -l doPseudoFits.C+\(\"/scratch/klawhorn/EffSysStore/etapt_0_0000.dat\",\"CB_MuSelEff/analysis/plots/etapt_0.root\",2,1,2,1,\"test\",0,0,\"/scratch/klawhorn/EWKAnaR12a/Efficiency/Zmm_MuSelEff/probe..root\"\)