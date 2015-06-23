#!/bin/bash

#cp ../TheoUnc/nnpdf30.txt .
#cp ../TheoUnc/mmht2014.txt .
#cp ../TheoUnc/ct14.txt .

sed -i '1,2d' nnpdf30.txt
sed -i '1,2d' mmht2014.txt
sed -i '1,2d' ct14.txt

paste -d " " headers_for_pdf_uncert.txt nnpdf30.txt mmht2014.txt ct14.txt | awk '{ print $1, $2, " & ", $4, " & ", $6, " & ", $8, "\\\\" }' > pdf_uncert_table.tex

