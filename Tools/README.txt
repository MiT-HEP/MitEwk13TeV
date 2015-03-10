MitEwk/Tools/README.txt - Jay Lawhorn 8/7/13

------| MERGE NTUPLES |------

MergeNtuples.C takes a single text file as input. The first line of the text file is your desired output file name. The names of the files you want to combine should follow each on their own line. Run like:

     root -l -q MergeNtuples.C+\(\"the_list.txt\"\)

------| MERGE JSON |------

combine_JSON.py is a python script that can "and", "or", or "sub" two json files containing luminosity sections. The syntax is as following:

    python combine_JSON.py -a file_one -b file_two -o output_file -r operator

where operator can be "and", "or" or "sub". 

MergeJSON.sh is a simple shell wrapper for combine_JSON.py that combines more than two JSON files at a time. The syntax is as following:

    source MergeJSON.sh input_list output_file

------| Get Luminosity and Pileup Info |------

Go to:

    https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiCalc

For latest version of lumiCalc. In summer 2013, I used

    python lumiCalc2.py overview - json.txt

To get luminosity information. It's going to take a while, just as a warning. Also, you want recorded luminosity, not delivered.

To get pile up information / Npu distribution for my data, I used:

    python pileupCalc.py -i json.txt --inputLumiJSON=/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions12/8TeV/PileUp/pileup_latest.txt --calcMode=true --minBiasXsec=69300 --maxPileupBin 50 --numPileupBins 500 output.root
