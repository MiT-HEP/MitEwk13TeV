#! /bin/bash

mkdir -p plots

python smoothSysUnfoldModel.py -c Zmm -v ZPt -s 1 -w 0.4
python smoothSysUnfoldModel.py -c Zmm -v PhiStar -s 1 -w 0.4
python smoothSysUnfoldModel.py -c Zmm -v ZRap -s 0 -w 0.4
python smoothSysUnfoldModel.py -c Zmm -v Lep1Pt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zmm -v Lep2Pt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zmm -v LepNegPt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zmm -v LepPosPt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zmm -v Lep1Eta -s 0 -w 0.3
python smoothSysUnfoldModel.py -c Zmm -v Lep2Eta -s 0 -w 0.3
