#!/bin/bash
#-----------------------------------------------------------------
# put this in POWHEG-BOX/W/ along with *-powheg.init files
# and execute from there
#-----------------------------------------------------------------

COM=13
BEAM_EN=6500
PDFSET=11000

mkdir wpm_${COM}
cp wpm-powheg.input wpm_${COM}/wpm-powheg.input
cd wpm_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" wpm-powheg.input
echo wpm | ../pwhg_main
cd ..

mkdir wmm_${COM}
cp wmm-powheg.input wmm_${COM}/wmm-powheg.input
cd wmm_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" wmm-powheg.input
echo wmm | ../pwhg_main
cd ..

mkdir wpe_${COM}
cp wpe-powheg.input wpe_${COM}/wpe-powheg.input
cd wpe_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" wpe-powheg.input
echo wpe | ../pwhg_main
cd ..

mkdir wme_${COM}
cp wme-powheg.input wme_${COM}/wme-powheg.input
cd wme_${COM}
sed -i "s/BEAM_EN/${BEAM_EN}/;s/PDFSET/${PDFSET}/" wme-powheg.input
echo wme | ../pwhg_main
cd ..

