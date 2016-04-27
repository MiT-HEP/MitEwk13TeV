#! /bin/bash

mkdir -p plots

python smoothSysUnfoldModel.py -c Zee -v ZPt -s 1 -w 0.4
python smoothSysUnfoldModel.py -c Zee -v PhiStar -s 1 -w 0.4
python smoothSysUnfoldModel.py -c Zee -v ZRap -s 0 -w 0.4
python smoothSysUnfoldModel.py -c Zee -v Lep1Pt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zee -v Lep2Pt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zee -v LepNegPt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zee -v LepPosPt -s 1 -w 0.2
python smoothSysUnfoldModel.py -c Zee -v Lep1Eta -s 0 -w 0.3
python smoothSysUnfoldModel.py -c Zee -v Lep2Eta -s 0 -w 0.3
