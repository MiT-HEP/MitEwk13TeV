import sys

def concatenate(listline,delimiter):
    stringline = ''
    for j in range(len(listline)):
        elt = listline[j]
        if j<len(listline)-1:
            stringline = stringline+elt+delimiter
        else:
            stringline = stringline+elt
    return stringline

efftype = sys.argv[1]
# File names to process. Always put the positive before the negative.
outputdir = '/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results-coarsebinning/'
infilenames = [outputdir+efftype+'/latex.txt']
# Output file
outfile = open(efftype+'_chargeIndependent_table.txt','w')

# Get number of eta bins (number of columns)
file_ = open(infilenames[0])
data = file_.readlines()
etabins = data[0]
needMultRows = False
netabins = etabins.count('eta')
if netabins>6:
	needMultRows = True
tablecolumns = 'cr'
for column in range(netabins-1):
	tablecolumns += '|c'
	# Muon trigger efficiency has too many eta bins to fit in just two rows
	if needMultRows and column>2:
		break
tablecolumns += '|cc'

# LaTex table setup
outfile.write('\\begin{table} \n')
outfile.write('\\begin{center} \n')
outfile.write('\\begin{tabular}{'+tablecolumns+'} \n')
outfile.write('\hline \n')

lines_firstrow = []
lines_secondrow = [] # For muon trigger efficiency

for filename in infilenames:
    # Open file
    infile = open(filename)
    datafile = infile.readlines()

    # Parse each line
    for linenum in range(len(datafile)):
        line = datafile[linenum]
        if 'eta' in line:
		if needMultRows:
			lineparsed = line.split('&')
			etabins_firstrow = concatenate(lineparsed[:(int)(netabins/2)+2],'&')+' \\\\'
			etabins_secondrow = '& &'+concatenate(lineparsed[(int)(netabins/2)+2:],'&')
		else:
			etabins_firstrow = line
			etabins_secondrow = None
        elif 'p_{T}' in line:
            	lineparsed = line.split('$')
            	ptbin = '$'+lineparsed[1]+'$&'#lineparsed[0]+'$'+lineparsed[1]+' (+)$'
		if needMultRows:
			effs_firstrow = ptbin+concatenate(lineparsed[2:netabins+2],'$')+'$ \\\\'
			effs_secondrow = ptbin+concatenate(lineparsed[netabins+2:],'$')
			lines_firstrow.append(effs_firstrow)
			lines_secondrow.append(effs_secondrow)
		else:
			effs_firstrow = ptbin+concatenate(lineparsed[2:],'$')
			lines_firstrow.append(effs_firstrow)

    infile.close()

outfile.write(etabins_firstrow)
outfile.write('\n \hline \hline \n')
for i in range(len(lines_firstrow)):
    outfile.write(lines_firstrow[i])
    outfile.write('\n \hline \n')
if needMultRows:
	outfile.write('\hline \n')
	outfile.write(etabins_secondrow)
	outfile.write('\n \hline \hline \n')
	for k in range(len(lines_secondrow)):
		outfile.write(lines_secondrow[k])
		outfile.write('\n \hline \n')

# Get table caption
efftype_caption_map = {'Zee_EleHLTEff':'Electron trigger efficiencies in MC.',
		       'DataZee_EleHLTEff':'Electron trigger efficiencies in data.',
		       'Zee_EleGsfSelEff':'Electron reconstruction, ID, and isolation efficiencies in MC.',
		       'DataZee_EleGsfSelEff':'Electron reconstruction, ID, and isolation efficiencies in data.',
		       'Zee_EleGsfEff':'Electron reconstruction efficiencies in MC.',
		       'DataZee_EleGsfEff':'Electron reconstruction efficiencies in data.',
		       'Zee_EleGsfEff_dRmatching':'Electron reconstruction (via deltaR matching) efficiencies in MC.',
		       'DataZee_EleGsfEff_dRmatching':'Electron reconstruction (via deltaR matching) efficiencies in data.',
		       'Zee_EleSCEff':'Electron supercluster efficiencies in MC.',
                       'Zmm_MuHLTEff':'Muon trigger efficiencies in MC.',
                       'DataZmm_MuHLTEff':'Muon trigger efficiencies in data.',
                       'Zmm_MuSelEff':'Muon ID and isolation efficiencies in MC.',
                       'DataZmm_MuSelEff':'Muon ID and isolation efficiencies in data.',
                       'Zmm_MuTrkEff':'Tracking efficiency in MC.',
                       'DataZmm_MuTrkEff':'Tracking efficiency in data.',
                       'Zmm_MuStaEff_iso':'Stand-alone muon efficiency in MC.',
                       'DataZmm_MuStaEff_iso':'Stand-alone muon efficiency in data.',
                       'DataZmm_MuStaEff_iso_extraTrkCuts':'Stand-alone muon efficiency in data.'}
caption = efftype_caption_map[efftype]

# LaTex table end
outfile.write('\n\end{tabular}\n')
outfile.write('\caption{\label{tab:'+efftype+'_table}'+caption+'}\n')
outfile.write('\end{center}\n')
outfile.write('\end{table}\n')
outfile.close()

