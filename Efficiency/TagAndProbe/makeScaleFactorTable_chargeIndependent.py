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
outputdir = '/afs/cern.ch/work/c/cmedlock/public/wz-efficiency-results/DataMC/'
print('***** USING EFFICIENCY RESULTS IN ',outputdir)
infilenames = [outputdir+efftype+'_scalefactors.txt']
# Output file
outfile = open(efftype+'_scalefactors_table.txt','w')

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
efftype_caption_map = {'eleHLTEff':'Electron trigger scale factors.',
                       'eleHLTEff_binsFromKevin':'Electron trigger scale factors.',
                       'eleHLTEff_run251244':'Electron trigger scale factors for run 251244.',
                       'eleHLTEff_run251251':'Electron trigger scale factors for run 251251.',
                       'eleHLTEff_run251252':'Electron trigger scale factors for run 251252.',
                       'eleHLTEff_run251561':'Electron trigger scale factors for run 251561.',
                       'eleHLTEff_run251562':'Electron trigger scale factors for run 251562.',
                       'eleHLTEff_run251643':'Electron trigger scale factors for run 251643.',
                       'eleHLTEff_run251721':'Electron trigger scale factors for run 251721.',
                       'eleHLTEff_run251883':'Electron trigger scale factors for run 251883.',
		       'eleGsfEff':'Electron reconstruction scale factors.',
                       'eleSCEff':'Electron supercluster scale factors.',
		       'eleGsfEff_dRmatching':'Electron reconstruction (via dR matching) scale factors.',
		       'eleSelEff':'Electron ID and isolation scale factors.',
		       'eleSelEff_binsFromKevin':'Electron ID and isolation scale factors.',
                       'eleSelEff_tightID':'Electron ID and isolation scale factors.',
                       'eleSelEff_run251244':'Electron ID and isolation scale factors for run 251244.',
                       'eleSelEff_run251251':'Electron ID and isolation scale factors for run 251251.',
                       'eleSelEff_run251252':'Electron ID and isolation scale factors for run 251252.',
                       'eleSelEff_run251561':'Electron ID and isolation scale factors for run 251561.',
                       'eleSelEff_run251562':'Electron ID and isolation scale factors for run 251562.',
                       'eleSelEff_run251643':'Electron ID and isolation scale factors for run 251643.',
                       'eleSelEff_run251721':'Electron ID and isolation scale factors for run 251721.',
                       'eleSelEff_run251883':'Electron ID and isolation scale factors for run 251883.',
		       'eleGsfSelEff':'Electron reconstruction, ID, and isolation scale factors.',
                       'muHLTEff':'Muon trigger scale factors.',
                       'muSelEff':'Muon ID and isolation scale factors.',
                       'muTrkEff':'Muon tracking scale factors.',
                       'muTrkEff_v2':'Muon tracking scale factors.',
                       'muSelStaEff':'Muon stand-alone, ID and isolation scale factors.',
                       'muSelStaEff_ratios':'Muon stand-alone, ID and isolation scale factors: relative difference between factorized and combined values, defined as ( SF(stand-alone)*SF(ID+Iso) - SF(stand-alone+ID+Iso) ) / SF(stand-alone+ID+Iso).',
                       'muStaEff':'Stand-alone muon efficiency scale factors.'}
caption = efftype_caption_map[efftype]

# LaTex table end
outfile.write('\n\end{tabular}\n')
outfile.write('\caption{\label{tab:'+efftype+'_table}'+caption+'}\n')
outfile.write('\end{center}\n')
outfile.write('\end{table}\n')
outfile.close()

