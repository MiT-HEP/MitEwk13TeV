#! /bin/bash

mkdir -p plots

python smoothSysUnfoldModel_ZPt.py
python smoothSysUnfoldModel_PhiStar.py
python smoothSysUnfoldModel_ZRap.py
python smoothSysUnfoldModel_Lep1Pt.py
python smoothSysUnfoldModel_Lep2Pt.py
python smoothSysUnfoldModel_LepNegPt.py
python smoothSysUnfoldModel_LepPosPt.py
python smoothSysUnfoldModel_Lep1Eta.py
python smoothSysUnfoldModel_Lep2Eta.py
