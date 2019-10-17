#! /bin/bash
## Take the plot and text file table output and attempt to make the latex files so I don't kill myself instead

# do the text files first....
INDIR="/afs/cern.ch/user/s/sabrandt/work/public/FilesSM2017GH/Efficiency/LowPU2017ID_13TeV/results"
OUTDIR="${INDIR}/latex"

SFX="_aMCxPythia"
COM="13"
# COM="5"


mkdir -p ${OUTDIR}/fig
mkdir -p ${OUTDIR}/tab
# mkdir -p ${OUTDIR}/plots/${COM}TeV

# my_array=(foo bar)
# for i in "${my_array[@]}"; do echo "$i"; done


LEPS=("Mu")
LEPWS=("muons")
Q=("Combined")
CAT_MU=("Sta")
PFX=("Zmm")
# LEPS=("Ele" "Mu")
# LEPWS=("electrons" "muons")
# Q=("Positive" "Negative")
# CAT_ELE=("GSFSel" "HLT")
# CAT_MU=("Sel" "Sta" "HLT")
# CAT_MU=("SIT" "HLT")
# PFX=("Zee" "Zmm")
SAMP=("Data" "MC")
TXTF="PtBins_scalefactors_niceTable.txt"
TABPFX="tab/tab.07.04.Eff"
FIGPFX="fig/fig.07.04.Eff"
## Loops to set up the directory paths 
CAP_ELE=("GSF electron identification and isolation" "Single electron trigger")
CAP_MU=("Muon selection" "Standalone muon identification" "Single muon trigger")
CAPQ=("positively" "negatively")
for ((i=0;i<${#LEPS[@]};++i)); do # BEGIN loop through lepton flavors
  flav=${LEPS[i]}
  echo ${flav}
  if [[ "${flav}" = "Ele" ]]; then
	pfx="Zee"
	CAT=("${CAT_ELE[@]}")
	CAP=("${CAP_ELE[@]}")
  elif [[ "${flav}" = "Mu" ]]; then
    pfx="Zmm"
	CAT=("${CAT_MU[@]}")
	CAP=("${CAP_MU[@]}")
  fi
  echo ${CAT}
  for ((j=0;j<${#CAT[@]};++j)); do
    cat=${CAT[j]}
    # Make a new .tex file which holds the table for the + and - tables for this bin
    TABLEFILE="${OUTDIR}/${TABPFX}.${flav}.${COM}TeV.${cat}.tex"
    echo "%%%% Tables for ${flav}${cat} Efficiency  %%%%%" > ${TABLEFILE}
    # Make a new .tex file that holds the figures for + and - for this thing
    PLOTSFILE="${OUTDIR}/${FIGPFX}.${flav}.${COM}TeV.${cat}.tex"
    echo "%%%% Figures for ${flav}${cat} Efficiency  %%%%%" > ${PLOTSFILE}
    
    echo ${TABLEFILE}
    for ((k=0;k<${#Q[@]};++k)); do # BEGIN loop through charges
      q=${Q[k]}
      PLOTDIR="figures/${COM}TeV/${pfx}${cat}${SFX}/${q}"
      PLOTDIR=${PLOTDIR,,}
      mkdir -p ${OUTDIR}/${PLOTDIR}
      echo ${OUTDIR}/${PLOTDIR}
      echo "Copying plots"
      cp ${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}/*.png ${OUTDIR}/${PLOTDIR}
      ls ${OUTDIR}/${PLOTDIR}
      
      #########  table shit  ###################
      echo "%% Efficiency table for ${pfx}${cat} ${q}" >> ${TABLEFILE}
      echo "directory name is  ${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}"
      echo "\begin{table}%[htbp]" >> ${TABLEFILE}
      echo "\begin{center}" >> ${TABLEFILE}
      echo "\scalebox{0.7}{" >> ${TABLEFILE}
      ## Add the file that has the formatting & data for the double-row table that doesn't look like shit
      cat "${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}/${TXTF}" >> ${TABLEFILE}
      echo "\end{center}" >> ${TABLEFILE}
      # echo "bbbbbbbbbbbbbb   ${CAP[@]}"
      echo "\caption{${CAP[j]} efficiency scale factors in (\$p_T\$, \$\eta\$) bins for ${CAPQ[k]} charged ${LEPWS[i]} in the ${COM} TeV samples.}" >> ${TABLEFILE}
      echo "\label{tab:Eff:${LEPWS[i]:0:2}:${COM}TeV:${cat}:${CAPQ[k]:0:3}}" >> ${TABLEFILE}
      echo "\end{table}" >> ${TABLEFILE}
      # cat ${TABLEFILE}
      ########################################################
      ###       figure shit     ############################
      # get a list of only the .png files in the directory
      
      echo "\begin{figure}" >> ${PLOTSFILE}
      PLOTS=( $( ls ${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}/ | grep ".png") )
      for ((m=0;m<${#PLOTS[@]};++m)); do
        echo ${PLOTS[m]}
        echo "\centering" >> ${PLOTSFILE}
        echo "\begin{subfigure}{.33\textwidth}" >> ${PLOTSFILE}
        echo "\centering" >> ${PLOTSFILE}
        PLOTNAME="ch07/figures/${COM}TeV/${pfx}${cat}${SFX}/${q}"
        PLOTNAME=${PLOTNAME,,} # lowercase that shit because Overleaf keeps ruining it randomly
        echo ${PLOTNAME}
        echo "\includegraphics[width=\linewidth]{${PLOTNAME}/${PLOTS[m]}}" >> ${PLOTSFILE}
        echo "%   \caption{1a}" >> ${PLOTSFILE}
        echo "%   \label{fig:Eff:${LEPWS[i]:0:2}:${COM}TeV:${cat}:${CAPQ[k]:0:3}" >> ${PLOTSFILE}
        echo "\end{subfigure}%" >> ${PLOTSFILE}
        echo $(( ${m} % 2 ))
        if [[ $(( ${m} % 3 )) = 0 ]] || [[ $(( ${m} % 3 )) = 1 ]]; then
          echo "hello"
          # echo "\hfill" >> ${PLOTSFILE}
        elif [[ $(( ${m} % 3 )) = 2 ]]; then
          echo "\\\\ " >> ${PLOTSFILE}
        fi
      done # END loop through pT bins
      echo "\caption{\$\eta\$ dependence of ${CAP[j]} efficiency scale factors, separated by \$p_T\$ bins, for ${CAPQ[k]} charged ${LEPWS[i]} in the ${COM} TeV samples.}" >> ${PLOTSFILE}
      echo "\end{figure}" >> ${PLOTSFILE}
      
      
	done # END loop through charges
  done # END loop through efficiency types
done # END loop through lepton flavors