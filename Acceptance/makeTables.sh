#! /bin/bash
## Take the plot and text file table output and attempt to make the latex files so I don't kill myself instead

# do the text files first....


SFX="_v1"

INDIR="."
OUTDIR="${INDIR}/latex${SFX}"
# COM="13TeV"
COM="5TeV"


mkdir -p ${OUTDIR}
# mkdir -p ${OUTDIR}/tab
# mkdir -p ${OUTDIR}/plots/${COM}

# my_array=(foo bar)
# for i in "${my_array[@]}"; do echo "$i"; done

TYPES=("${SFX}" "_dressed_${SFX}")
LEP_OPTS=("Wep" "Wem" "We" "Zee" "Wmp" "Wmm" "Wm" "Zmm")
LABELS=("$W\rightarrow e^+\nu$" "$W\rightarrow e^-\nu$" "$W\rightarrow e\nu$" "$Z\rightarrow ee$" "$W\rightarrow \mu^+\nu$" "$W\rightarrow \mu^-\nu$" "$W\rightarrow \mu\nu$" "$Z\rightarrow \mu\mu$" )

TXTF="gen.txt"
TABPFX="tab.09.01.Acc.Gen"
TABLEFILE="${OUTDIR}/${TABPFX}.Val.${COM}.tex"
echo "%%%% Tables for the ${COM}"  > ${TABLEFILE}
## Loops to set up the directory paths 
      echo "%% Efficiency table for ${pfx}${cat} ${q}" >> ${TABLEFILE}
      # echo "directory name is  ${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}"
      echo "\begin{table}%[htbp]" >> ${TABLEFILE}
      echo "\begin{center}" >> ${TABLEFILE}
      echo "\scalebox{0.7}{" >> ${TABLEFILE}
for ((i=0;i<${#LEP_OPTS[@]};++i)); do # BEGIN loop through the dressed/undressed types
  lep=${LEP_OPTS[i]}
  echo ${dress}
  for ((j=0;j<${#TYPES[@]};++j)); do # BEGIN loop through the options
    # cat=${CAT[j]}
    # Make a new .tex file which holds the table for the + and - tables for this bin
    echo ${TABLEFILE}

      ## Add the file that has the formatting & data for the double-row table that doesn't look like shit
      cat "${INDIR}/Plots/${pfx}${cat}${SFX}/pt/${q}/${TXTF}" >> ${TABLEFILE}

      # cat ${TABLEFILE}
      ########################################################
      ###       figure shit     ############################
      # get a list of only the .png files in the directory
      
      
	# done # END loop through charges
  done # END loop through Boson/lepton configurations
  
  
done # END loop through dressed/undressed
echo "\end{center}" >> ${TABLEFILE}
# echo "bbbbbbbbbbbbbb   ${CAP[@]}"
echo "\caption{Acceptance values for the post-FSR and dressed leptons for ${COM}.}" >> ${TABLEFILE}
echo "\label{tab:Acc:Gen:Val:${COM}}}" >> ${TABLEFILE}
echo "\end{table}" >> ${TABLEFILE}