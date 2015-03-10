#!/usr/bin/env python

import sys,ConfigParser,os,string,commands,time
from optparse import OptionParser
import json
from FWCore.PythonUtilities.LumiList import LumiList


def main():
    
    parser=OptionParser()
    # required parameters
    parser.add_option('-a','--json1',dest='JSONFILE1',default='none',help='first JSON file')
    parser.add_option('-b','--json2',dest='JSONFILE2',default='none',help='second JSON file')
    parser.add_option('-o','--output',dest='OUTPUT',default='none',help='output JSON file')
    parser.add_option('-r','--operation',dest='OPERATION',default='none',help='and,or,sub')
    
    (options, args) = parser.parse_args()

    if options.JSONFILE1=='none':
        print "Please provide a json1 file"
        sys.exit()
    if options.JSONFILE2=='none':
        print "Please provide a json2 file"
        sys.exit()
    if options.OUTPUT=='none':
        print "Please provide a json2 file"
        sys.exit()
    if options.OPERATION!='and' and options.OPERATION!='or' and options.OPERATION!='sub':
        print "Please provide an operation (and,or,sub)"
        sys.exit()

    J1=options.JSONFILE1
    J2=options.JSONFILE2
    OUTFILE=options.OUTPUT
    op=options.OPERATION

    J1List = LumiList(filename = J1)  # Read in first JSON file
    J2List = LumiList(filename = J2)  # Read in first JSON file

#    print J1List
    OutList={}
    if op=='and':
        OutList =J1List & J2List
    if op=='or':
        OutList =J1List | J2List
    if op=='sub':
        OutList =J1List - J2List
            
#    print OutList.getCompactList()

# Write out the output list

    lumiout = open(OUTFILE, 'w')
    json.dump(OutList.getCompactList(), lumiout,sort_keys=True)
    lumiout.close()


if __name__ == '__main__':
    main()
