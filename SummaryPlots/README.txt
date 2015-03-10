MitEwk/SummaryPlots/README.txt - Jay Lawhorn 8/8/13

------| INTRODUCTION |------

This module takes the values measured in all the other modules and calculates the cross section of the different channels. It also compares these values
to the theoretical predictions.

------| RUN |------

* root -l -q xsec.C+\(\"input.txt\", \"output_dir\"\)

------| MODIFICATION |------

* in xsec.C:

    theoretical predictions for the various cross sections are defined

* in input.txt:

  yield_err and acc_x_sf_err are absolute

  all other uncertainties are relative

  luminosity uncertainty should be obtained from the dedicated group for that