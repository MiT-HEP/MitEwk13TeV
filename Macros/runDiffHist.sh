#! bin/bash

source roofitsetup.sh

printf 'making plots\n'
#root -l -q getDiffHist.C+\(\"../Zee/mlfit_zee.root\",\"Zee\"\) && mv *.png /home/cmedlock/public_html/Combine/Zee/plots
root -l -q getDiffHist.C+\(\"../Zmm/mlfit_zmm.root\",\"Zmm\"\) && mv *.png /home/cmedlock/public_html/Combine/Zmm/plots

#root -l -q getDiffHist.C+\(\"../Wenu/mlfit_wenu.root\",\"Wenu\"\)
#root -l -q getDiffHist.C+\(\"../Wenu/mlfit_wenu_p.root\",\"Wenu_p\"\)
#root -l -q getDiffHist.C+\(\"../Wenu/mlfit_wenu_m.root\",\"Wenu_m\"\)
#mv *.png /home/cmedlock/public_html/Combine/Wenu/plots

#root -l -q getDiffHist.C+\(\"../Wmunu/mlfit_wmunu.root\",\"Wmunu\"\)
#root -l -q getDiffHist.C+\(\"../Wmunu/mlfit_wmunu_p.root\",\"Wmunu_p\"\)
#root -l -q getDiffHist.C+\(\"../Wmunu/mlfit_wmunu_m.root\",\"Wmunu_m\"\)
#mv *.png /home/cmedlock/public_html/Combine/Wmunu/plots

printf 'cleaning up\n'
rm *_C.d *_C.so *~
