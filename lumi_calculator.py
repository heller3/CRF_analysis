#!/usr/bin/env python
import os, sys, re
import os
from array import array

## This function takes RunList.txt and FillReport.txt and returns a list of integrated luminosities for each run
outFile = open("LumiList.txt","w")

with open("RunList.txt") as runlist:
    for runnum in runlist:
        print "Run ", str(int(runnum))
        totallumi=0.
        totaltime=0.
        with open("FillReport.txt") as fills:
            for fill in fills:
                if "Fill" in fill.split("\t")[0]: continue
                if int(fill.split("\t")[24].split()[0]) < int(runnum):
                    print "\tFill ",fill.split("\t")[0], ", lumi ",fill.split("\t")[7],", duration ",fill.split("\t")[2]
                    print "\tRun list is " ,fill.split("\t")[24].split()
                    #heavy ion runs are measured in inverse microbarns, not pico
                    if "Ion" not in fill.split("\t")[14]:
                        totallumi = totallumi+float(fill.split("\t")[7])
                    else:
                        totallumi = totallumi+(float(fill.split("\t")[7])/1000000.)
                    totaltime = totaltime+ 60*int(fill.split("\t")[2].split()[0])+int(fill.split("\t")[2].split()[2])

        totaltime=totaltime/60.
        print "Total Lumi: ",totallumi/1000.," fb-1"
        print "Total dose time: ",totaltime," hours"
        print "Average dose rate: ",totallumi/totaltime," pb-1 per hour"

        outFile.write("%s %s %s" %(str(int(runnum)),totallumi/1000., str(totallumi/totaltime)+"\n"))
        
outFile.close()
