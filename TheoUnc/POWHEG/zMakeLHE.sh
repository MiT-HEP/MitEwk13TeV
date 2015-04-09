#!/bin/bash
#-----------------------------------------------------------------
# put this in POWHEG-BOX/Z/ along with *-powheg.init files
# and execute from there
#-----------------------------------------------------------------

COM=13
BEAM_EN=6500
PDFSET=11000

mkdir zmm_${COM}
cp zmm-powheg.input zmm_${COM}/zmm-powheg.input
cd zmm_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" zmm-powheg.input
echo zmm | ../pwhg_main
cd ..

mkdir zee_${COM}
cp zee-powheg.input zee_${COM}/zee-powheg.input
cd zee_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" zee-powheg.input
echo zee | ../pwhg_main
cd ..

