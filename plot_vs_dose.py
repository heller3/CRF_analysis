#!/usr/bin/env python
import os, sys, re
import ROOT
import os
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetLabelFont(42,"xyz")
ROOT.gStyle.SetLabelSize(0.05,"xyz")
#ROOT.gStyle.SetTitleFont(42)
ROOT.gStyle.SetTitleFont(42,"xyz")
ROOT.gStyle.SetTitleFont(42,"t")
#ROOT.gStyle.SetTitleSize(0.05)
ROOT.gStyle.SetTitleSize(0.06,"xyz")
ROOT.gStyle.SetTitleSize(0.06,"t")
ROOT.gStyle.SetPadBottomMargin(0.14)
ROOT.gStyle.SetPadLeftMargin(0.14)
ROOT.gStyle.SetTitleOffset(1,'y')
#ROOT.gStyle.SetLegendTextSize(0.05)
ROOT.gStyle.SetGridStyle(3)
ROOT.gStyle.SetGridColor(13)

one = ROOT.TColor(2001,0.906,0.153,0.094)
two = ROOT.TColor(2002,0.906,0.533,0.094)
three = ROOT.TColor(2003,0.086,0.404,0.576)
four =ROOT.TColor(2004,0.071,0.694,0.18)
five =ROOT.TColor(2005,0.388,0.098,0.608)
six=ROOT.TColor(2006,0.906,0.878,0.094)
colors = [1,2001,2002,2003,2004,2005,2006,2,3,4,6,7,5,1,8,9,29,38,46]


#######These are default parameters than can be overwritten at command line

#usage
##./plot_vs_time.py <number> <run type> <axis type>

#number: indicates the set of channels used (defined below)
#run type: "ped" switches to pedestal mode; default is laser mode
#axis type: "Mrad" switches puts dose on x-axis; default is fb-1

#to distinguish laser mode or pedestal drift mode
laser = True

#to distinguish fb-1 or Mrad on x-axis
ifb = True

#to enable shunt correction
fixshunt=False

## to plot damage or recovery
recovery=False


shuntval = 3.87

#to apply a cut to runs considered (not ironed out yet)
applycut=False


if len(sys.argv) < 4:
    sys.exit('Please provide selection number, laser or ped mode, and xaxis mode') 

else:
    variation = int(sys.argv[1])
    if "ped" in sys.argv[2]: laser = False
    if "Mrad" in sys.argv[3]: ifb = False
    if "recovery" in sys.argv[3]: recovery=True

if len(sys.argv) > 4:
    if "unshunt" in sys.argv[4]: fixshunt=True
if len(sys.argv) > 5:
    if "cut" in sys.argv[5]: applycut=True        

#####################        


###### Run and luminosity lists

#histname defines which energy definition is used to define the mean
#Lasers use peak finding 4TS energy (h1_energy_ieta)
#Pedestal drift measurement uses constant 8TS energy (h1_energy_const_ieta)


#lasers
if laser:
    histname = "h1_energy_laserconst_ieta"
    if not recovery:
        runlist = [280662,281740,282046,282368,282369,282779,282933,283403,283512,283660,283803, 283842, 283973, 284127]
        lumi =    [0.    ,1.1   ,2.31   ,2.70  ,2.72,  3.29, 4.01  ,5.90  ,6.8  ,7.0   ,7.24  , 7.66  , 9.53  , 10.0 ]

    else:
        runlist = [280662, 284127, 284189, 284342, 284493,285227,285356,  286185, 286395 ]
        lumi =    [-100, 2     , 3       ,  6  , 8     ,18  ,    20 ,  35    , 37     ]
        # 284248, 284314,
         #5     ,   6   ,
   # runlist = [280662,281740,282046,282368,282369,282779,282933,283403,283512,283660,283803, 283842, 283918, 283973, 284127,284284,284342,284493,285356,286185]
    #lumi =    [0.    ,1.1   ,2.31   ,2.70  ,2.72,  3.29, 4.01  ,5.90  ,6.8  ,7.0   ,7.24  , 7.66  , 8.81  , 9.53  , 10.0    ,11   , 12    , 13   ,14    ,15   ]
    
    #lumi =    [0.    ,1.1   ,2.1   ,2.49  ,2.51  ,  3.  ,3.73  ,5.61  ,6.52  ,6.73  ,6.96  , 7.37  , 8.53  , 9.25] original lumi calculation


    ## currently exclude 283918 since it does not have a matching pedestal measurement
    
#pedestals
else:
    histname = "h1_energy_const_ieta"
    if not recovery:
        runlist = [280646,281734,282040,282314,282772,282927,283321,283397,283522, 283662,283807, 283846, 283977, 284123]
        lumi =    [0.    ,1.1   ,2.31   ,2.71,  3.29, 4.01  , 5.44 ,5.90  ,6.8  ,7.0   ,7.24  , 7.66  , 9.53  , 10.0 ]
        shunted = [0     ,0     ,1     ,1    ,1  ,     1    ,1     ,1     , 1      ,1      ,0     , 0    , 0,     0]

    ## recovery set
    else:
        runlist = [284123, 284185, 284249, 284306, 284337, 284487,285222,285347,  286184, 286394 ]
        lumi =    [2     , 3     , 5     ,   6   ,  6.3  , 8     ,18  ,    20 ,  35    , 37     ] ##285841 at 30 days, 286457 is 38 284876 is 13 days
    
    #lumi =    [0.    ,1.1   ,2.1   ,2.5 ,  3.,    3.73  ,5.16  ,5.61  , 6.52  , 6.73  ,6.96  , 7.37  , 9.25] original lumi calculation
    # indicates whether shunt was activated for run or not


legtype = "dist"    
outname = histname.strip('_ieta').replace("h1_","")
if fixshunt: outname = outname+"_shunt_corr"
if applycut and laser: outname = outname+"_cut"
        
if len(runlist) > 10:
    outname = outname+"_v14_"


#### Define which set of channels to compare ("variation"), their positions and a name for the set
        
if variation==0:
    name="RM1_SCSN81-S_close"
    title="RM1 SCSN81-S, close to beam"
    channels =  [(2,4),(2,1),(5,0),(7,5),(3,0)]
    positions = ["Ref","2-2","2-1","1-2","1-1"]


if variation==1:
    name="RM1_SCSN81-S_far"
    title="RM1 SCSN81-S, far from beam"
   # channels =  [(2,4),(7,2),(8,5),(3,2),(9,0)]
   # positions = ["Ref","3-2","3-2","3-1","3-1"]
    channels =  [(2,4),(7,2),(8,5),(3,2)]
    positions = ["Ref","3-2","3-2","3-1"]
    
    #fibers = [2,7,8,3]
    #channels=[4,2,5,0]

if variation==2:
    name="RM2_SCSN81-S_close"
    title="RM2 SCSN81-S, close to beam"
    channels=   [(16,5),(15,1),(20,4),(20,5),(20,1)]
    positions = ["Ref" ,"2-2" ,"2-2" ,"2-2" ,"1-1"]
    #channels=   [(16,5),(15,1),(20,4),(20,5),(17,4),(21,4),(20,1)]
    #positions = ["Ref" ,"2-2" ,"2-2" ,"2-2" ,"2-1" ,"1-2" ,"1-1"]

if variation==3:
    name="RM2_SCSN81-S_far"
    title="RM2 SCSN81-S, far from beam"
    #channels=   [(16,5),(14,3),(16,1),(14,4),(15,4),(17,2)]
    #positions = ["Ref" ,"3-2" ,"3-2" ,"3-1" ,"3-1" ,"3-1" ]
    channels=   [(16,5),(14,3),(16,1),(14,4),(15,4)]
    positions = ["Ref" ,"3-2" ,"3-2" ,"3-1" ,"3-1" ]    


    
if variation==4:
    name="ref_comp"
    title="Reference tile comparison"
    channels =  [(16,5),(2,4)]
    positions = ["Ref. 2","Ref. 1"]


if variation==5:
    name="RM1_SCSN81-S_close_ref2"
    title="RM1 SCSN81-S, close to beam"
    channels =  [(16,5),(2,1),(5,0),(7,5),(3,0)]
    positions = ["Ref","2-2","2-1","1-2","1-1"]


if variation==6:
    name="RM1_SCSN81-S_far_ref2"
    title="RM1 SCSN81-S, far from beam"
   # channels =  [(2,4),(7,2),(8,5),(3,2),(9,0)]
   # positions = ["Ref","3-2","3-2","3-1","3-1"]
    channels =  [(16,5),(7,2),(8,5),(3,2)]
    positions = ["Ref","3-2","3-2","3-1"]

if variation==7:
    name="empty_channels"
    title="Empty Channels"
    channels = [(4,5),(6,5),(16,4),(19,3)]
    positions = ["empty","empty","empty","empty"]


if variation==8:
    name="EJ260S"
    title = "EJ260-S channels"
    channels= [(16,5),(3,1),(8,2),(5,5)]
    positions=["Ref","2-1","1-1","1-1" ]


    
if variation==9:
    name="EJ200_scintX"
    title = "Scintillator X and EJ200-S channels"
    channels= [(16,5),(19,1),(15,5),(18,5)]
    positions=["Ref","1-2","3-2","1-1" ]

if variation==10:
    name="EJ260_SCSN"
    title = "18-19 cm from beamline"
    channels= [(16,5),(3,1),(5,0),(2,1),(20,5)]
    positions=["Ref","2-1","2-1","2-2","2-2"]
    legtype = "mat"
    graphlabel = ["Ref","EJ260-S","SCSN81-S","SCSN81-S","SCSN81-S"]
#    fibers = [16,14,16,20]
#    channels=[5,3,1,1]


if variation==11:
    name="verycloseSCSN81"
    title = "11-14 cm from beamline"
    channels= [(16,5),(7,5),(20,1),(3,0)]
    positions=["Ref","1-2","1-1","1-1"]
    graphlabel = ["Ref","SCSN81-S","SCSN81-S","SCSN81-S"]
    legtype = "mat"

if variation==12:
    name="mediumSCSN81"
    title = "18-19 cm from beamline"
    channels= [(16,5),(5,0),(2,1),(20,5),(20,4)]
    positions=["Ref","2-1","2-2","2-2","2-2"]
    graphlabel = ["Ref","SCSN81-S","SCSN81-S","SCSN81-S","SCSN81-S"]
    legtype = "mat"        

if variation==13:
    name="farSCSN81"
    title = "24-27 cm from beamline"
    channels= [(16,5),(14,4),(15,4),(3,2),(8,5)]
    positions=["Ref","3-1","3-1","3-1","3-2"]
    graphlabel = ["Ref","SCSN81-S","SCSN81-S","SCSN81-S","SCSN81-S"]
    legtype = "mat" 

if variation==14:
    name="bright"
    title="SCSN81-S channels"
    channels= [(16,5),(14,3),(14,4),(2,1),(20, 4),(20,1)]
    positions=["Ref","3-2" ,"3-1","2-2" ,"2-2"  ,"1-1"]

if variation==15:
    name="EJ260_SCSN2"
    title = "EJ260"
    channels= [(16,5),(20,1),(3,0),(20,4),(14,4),(3,1)]
    positions=["Ref","1-1","1-1","2-2","3-1","2-2"]
    legtype = "mat"
    graphlabel = ["Ref","SCSN81-S","SCSN81-S","SCSN81-S","SCSN81-S","EJ260-S"]


if variation==16:
    name="EJ200_LS"
    title = "Other"
    channels= [(16,5),(19,1),(16,2),(18,5),(15,5)]
    positions=["Ref","1-2","3-3","1-1","3-2"]
    legtype = "mat"
    graphlabel = ["Ref","Scintillator X","LS","EJ200-S","EJ200-S"]    

if variation==17:
    name="petpen"
    title = "PET/PEN/PTP"
    channels= [(16,5),(21,3),(21,5),(21,1),(21,2)]
    positions=["Ref","1-2","1-2","1-1","1-1"]
    legtype = "mat"
    graphlabel = ["Ref","PEN","PEN","PTP","PTP"]    

if variation==18:
    name="SCSN81F"
    title = "SCSN81-F"
    channels= [(16,5),(4,3),(15,2),(15,3),(15,0)]
    positions=["Ref","1-3","2-3","2-3","3-3"]
    #legtype = "mat"
    graphlabel = ["Ref","PEN","PEN","PTP","PTP"]    
    
if variation==19:
    name="bright_variety"
    title = "Other materials"
    channels= [(16,5),(4,3),(15,2),(3,1),(15,5)]
    positions=["Ref","1-3","2-3",  "2-2","3-2"]
    legtype = "mat"
    graphlabel = ["Ref","SCSN81-F","SCSN81-F","EJ260-S","EJ200-S"]    

if variation==20:
    name="scsnvar"
    title="SCSN81-S channels"
    channels= [(16,5),(14,4),(15,4),(20,4),(20,1),(3,0)]
    positions=["Ref","3-1" ,"3-1","2-2" ,"1-1" ,"1-1"]    
    graphlabel = ["","","","","",""]

if variation==21:
    name="othervar"
    title = "Other materials"
    channels= [(16,5),(4,3),(15,2),(3,1),(18,5)]
    positions=["Ref" ,"1-3","2-3" ,"2-2","1-1"]
    legtype = "mat"
    graphlabel = ["Ref","SCSN81-F","SCSN81-F","EJ260-S","EJ200-S"]   

### Map from position index to radius from beamline
    
dist_map = {"3-3":"27.2 cm","3-2":"25.9 cm","3-1":"24.6 cm",
            "2-3":"20.8 cm","2-2":"19.5 cm","2-1":"18.2 cm",
            "1-3":"14.4 cm","1-2":"13.1 cm","1-1":"11.8 cm",
            "empty":"-" }


### Map from position index to dose [Mrad per 4 fb-1]    
dose_map = {"3-3":0.57,"3-2":0.61,"3-1":0.7,
            "2-3":0.82,"2-2":0.85,"2-1":1.,
            "1-3":1.95,"1-2":2.21,"1-1":2.41,
            "empty":1.,"Ref":2.41 }

    




######## Measure QIE pedestal from run 280646, if needed to perform shunt correction #####
qie_currents = []
if fixshunt:
    pedFileName = "hists/hists_280646.root"
    pedFile = ROOT.TFile(pedFileName,"READ")
    for fiber,channel in channels:
        print "Getting hist "+histname+str(fiber-2)+"_iphi"+str(channel) 
        hist = pedFile.Get("h1_energy_const_ieta"+str(fiber-2)+"_iphi"+str(channel))
        print "Mean of raw hist: ",hist.GetMean()
        ##Find location of 0 SPE peak, and restrict range to 2 before and 2 after peak
        maxbin = hist.GetMaximumBin()
        hist.GetXaxis().SetRange(maxbin-2,maxbin+2)
        qie_currents.append(hist.GetMean()/8.) ### store in fC/25 ns
        print "Mean of 0 SPE peak: ", hist.GetMean()

        
leg_y_start = 0.63
if len(channels) <=2: leg_y_start = 0.7

v_mean = []
v_error = []
v_rms = []
v_uf = []

#### Iterate through runs and gather means
for run in runlist:
    print "run",str(run)
    inFileName = "hists/hists_"+str(run)+".root"
    inputFile = ROOT.TFile(inFileName,"READ")
    means = []
    errors = []
    rmss = []
    ufs = []
    for fiber,channel in channels:
        print str(fiber),str(channel)
        print "Getting hist "+histname+str(fiber-2)+"_iphi"+str(channel) 
        hist = inputFile.Get(histname+str(fiber-2)+"_iphi"+str(channel))
      
        mean = hist.GetMean()
        error = hist.GetMeanError()
        rms = hist.GetRMS()
        if hist.GetEntries()>0: uf = float(hist.GetBinContent(0))/hist.GetEntries()
        else: uf = 0.
        if "_const" in histname: ### this indicates pedestal run, divide by 8 to convert to single TS value
            mean = 0.125*mean
            error = 0.125*error
            rms = 0.125*rms
        print "mean is "+str(mean)
        #hist_ped = inputFile.Get("h1_fC_ieta"+str(fiber-2)+"_iphi"+str(channel))
       # if hist.GetBinContent(0)/hist.GetEntries() < 0.05:        
        means.append(mean)
        errors.append(error)
        rmss.append(rms)
        ufs.append(uf)
        # try storing RMS /underflow as well
        
             
    v_mean.append(means)
    v_error.append(errors)
    v_rms.append(rmss)
    v_uf.append(ufs)


#### Loop to produce graphs of mean fC (not normalized)
max = 1
graphs = []
for i,fiber in enumerate(channels):
    print str(i), str(fiber)
    x = []
    ex = []
    y = []
    ey = []
    for j in range(len(lumi)):
        #if fiber[0]==16 and fiber[1]==5 and j==0: continue
        if applycut and v_uf[j][i] > 0.05: continue
        thismean = v_mean[j][i]
        thiserr = v_error[j][i]
        if ((not laser) and fixshunt):
            if shunted[j] == 1:
              # thismean = shuntval*(thismean - qie_currents[i])+qie_currents[i]
              thismean = shuntval*(thismean - 18.75)+18.75
              thiserr = thiserr*shuntval

        y.append(thismean)    
        if thismean > max: max=v_mean[j][i]
        ey.append(thiserr)
        if ifb:
            x.append(lumi[j])
            ex.append(0.05)
        else:
            x.append(0.25*lumi[j]*dose_map[positions[i]])
            ex.append(0.05)    

    graphs.append(ROOT.TGraphErrors(len(x),array("d",x),array("d",y),array("d",ex),array("d",ey)))


c = ROOT.TCanvas()
c.SetGridy()
if laser:
    leg = ROOT.TLegend(0.47,leg_y_start,0.89,0.88)
elif recovery:
    leg = ROOT.TLegend(0.64,leg_y_start,0.89,0.88)
else:
    leg = ROOT.TLegend(0.14,0.6,0.37,0.88)

leg.SetMargin(0.15)
for i,graph in enumerate(graphs):
    #if not ifb:
        #print "triggered"
        #graph.GetXaxis().SetLimits(-0.6,6)
    graph.SetMinimum(0)
   # if "const" in histname:
    if laser: max = 2* max
    #elif recovery: max=350
    else: max = 300

    if laser and recovery: graph.GetXaxis().SetLimits(0,40)
    graph.SetMaximum(max)
    graph.SetLineColor(colors[i])
    graph.SetMarkerColor(colors[i])
    graph.SetMarkerSize(1)
    graph.SetMarkerStyle(20)
    if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Average charge [fC]")
    else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Average charge [fC]")
    if recovery: graph.SetTitle(title+"; Time elapsed since end of pp collisions [days]; Average charge [fC]")
    #leg.AddEntry(graph, "Chan. "+str(channels[i])+", position: "+positions[i],"EP")
    #if "const" in histname:
     #   leg.AddEntry(graph, "Chan. "+str(channels[i]),"EP")
    #else:     
     #   if "Ref" in positions[i]:  leg.AddEntry(graph, "Chan. "+str(channels[i])+", Reference tile","EP")
      #  else: leg.AddEntry(graph, "Chan. "+str(channels[i])+", Beam dist: "+dist_map[positions[i]],"EP")
    if laser:
        if "dist" in legtype:
            if "Ref" in positions[i]:  leg.AddEntry(graph, "Reference tile "+str(channels[i]),"EP")
            else: leg.AddEntry(graph, "R = "+dist_map[positions[i]]+"  "+str(channels[i]),"EP")
        else:
            if i>0: leg.AddEntry(graph, graphlabel[i]+" "+str(channels[i]),"EP")
    else:
        leg.AddEntry(graph, "Channel "+str(channels[i]),"EP")
        #else: leg.AddEntry(graph, "R = "+dist_map[positions[i]]+"  "+str(channels[i]),"EP")
            
    #graph.GetXaxis().SetLimits(-0.6,6)
    if i==0: graph.Draw("AELPZ")
    else: graph.Draw("ELPZ same")

leg.Draw()
if recovery:
     c.Print("doseplots/recovery_"+outname+"_"+name+"_vs_days.pdf")
else:           
    if ifb: c.Print("doseplots/"+outname+"_"+name+"_vs_lumi.pdf")
    else: c.Print("doseplots/"+outname+"_"+name+"_vs_dose.pdf")

    
if not laser: exit()

    
c.SetLogy();
for i,graph in enumerate(graphs):
    
    graph.SetMinimum(9)
    graph.SetMaximum(80*max)
   # if "const" in histname:
    #max = 2* max
    if i==0: graph.Draw("AELPZ")
    else: graph.Draw("ELPZ same")
leg.Draw()       

if recovery:
     c.Print("doseplots/log_recovery_"+outname+"_"+name+"_vs_days.pdf")
else:
    if ifb: c.Print("doseplots/log_"+outname+"_"+name+"_vs_lumi.pdf")
    else: c.Print("doseplots/log_"+outname+"_"+name+"_vs_dose.pdf")
    
        
c.Close()



#### Loop to form normalized graphs
normidx=[0]
for idx in normidx:
   
    max = 2
    #if "RM1" in name: max = 4
    graphs_norm = []
    c = ROOT.TCanvas()
    c.SetGridy();
    c.SetLogy();
    for i,fiber in enumerate(channels):
        print str(i), str(fiber)
        x = []
        ex = []
        y = []
        ey = []
        for j in range(len(lumi)):
            #if fibers[idx]==16 and j==0: continue
            if applycut and v_uf[j][i] > 0.05: continue
            y.append(v_mean[j][i]/v_mean[j][idx])
            if v_mean[j][i]/v_mean[j][idx] > max: max = v_mean[j][i]/v_mean[j][idx]
            ey.append(pow(pow(v_error[j][i],2)+pow(v_error[j][idx],2),0.5)/v_mean[j][idx])
            if ifb:
                x.append(lumi[j])
                ex.append(0.05)
            else:
                x.append(0.25*lumi[j]*dose_map[positions[i]])
                ex.append(0.05)    

        graphs_norm.append(ROOT.TGraphErrors(len(x),array("d",x),array("d",y),array("d",ex),array("d",ey)))

    leg = ROOT.TLegend(0.47,leg_y_start,0.89,0.89)
    leg.SetMargin(0.15)
    for i,graph in enumerate(graphs_norm):
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6.6)
        if laser and recovery: graph.GetXaxis().SetLimits(0,40)
        graph.SetMinimum(0.02)
        graph.SetMaximum(20)
        graph.SetLineColor(colors[i])
        graph.SetMarkerColor(colors[i])
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Average charge / reference charge")
        else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Average charge / reference charge")     #   graph.SetTitleFont(42);
        if recovery: graph.SetTitle(title+"; Time elapsed since end of pp collisions [days]; Average charge / reference charge")
        #leg.AddEntry(graph, "Chan. "+str(channels[i])+", position: "+positions[i],"EP")
        if "dist" in legtype:
            if "Ref" not in positions[i]: leg.AddEntry(graph, "R = "+dist_map[positions[i]]+"  "+str(channels[i]),"EP")
        else:
            if i>0: leg.AddEntry(graph, graphlabel[i]+" "+str(channels[i]),"EP")
                
        if i==1: graph.Draw("AELPZ")
        elif i>1: graph.Draw("ELPZ same")

    leg.Draw()
    if recovery:
        c.Print("doseplots/log_recovery_"+outname+"_"+name+"_vs_days_norm.pdf")
    else:
        if ifb: c.Print("doseplots/log_"+outname+"_"+name+"_vs_lumi_norm.pdf")
        else: c.Print("doseplots/log_"+outname+"_"+name+"_vs_dose_norm.pdf")
    c.Close()

    


### Loop again to form graph normalized, relative to initial value    
normidx=[0]
for idx in normidx:
    max = 2
    
    graphs_norm = []
    c = ROOT.TCanvas()
    c.SetGridy();
    for i,fiber in enumerate(channels):
        print str(i), str(fiber)
        x = []
        ex = []
        y = []
        ey = []
        y_first=1
        for j in range(len(lumi)):
           
            yval = v_mean[j][i]/v_mean[j][idx]
            if j==0: y_first = yval
            yval = yval/y_first
            if applycut and v_uf[j][i] > 0.05: continue
           # if v_mean[j][i] < 75:
            #    continue
            if recovery and j==0: continue
            y.append(yval)
            if yval > max: max = yval
            ey.append(pow(pow(v_error[j][i],2)+pow(v_error[j][idx],2),0.5)/v_mean[j][idx])
            if ifb:
                x.append(lumi[j])
                ex.append(0.05)
            else:
                x.append(0.25*lumi[j]*dose_map[positions[i]])
                ex.append(0.05)    


        graphs_norm.append(ROOT.TGraphErrors(len(x),array("d",x),array("d",y),array("d",ex),array("d",ey)))
    if "v9" not in outname:
       # if not recovery:
        leg = ROOT.TLegend(0.47,leg_y_start,0.89,0.89)
        #else:
        #leg = ROOT.TLegend(0.53,0.2,0.89,0.45)
    else:
        leg = ROOT.TLegend(0.65,0.65,0.89,0.89)
    leg.SetMargin(0.15)
    for i,graph in enumerate(graphs_norm):
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6.6)
        if laser and recovery: graph.GetXaxis().SetLimits(0,40)
        graph.SetMinimum(0)
        graph.SetMaximum(1.15*max)
        graph.SetLineColor(colors[i])
        graph.SetMarkerColor(colors[i])
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Remaining relative brightness")
        else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Remaining relative brightness")
        if recovery: graph.SetTitle(title+"; Time elapsed since end of pp collisions [days]; Remaining relative brightness")

        #   graph.SetTitleFont(42);
        #leg.AddEntry(graph, "Chan. "+str(channels[i])+", position: "+positions[i],"EP")
        if "dist" in legtype:
            if "Ref" not in positions[i]: leg.AddEntry(graph, "R = "+dist_map[positions[i]]+"  "+str(channels[i]),"EP")
        else:
            if i>0: leg.AddEntry(graph, graphlabel[i]+" "+str(channels[i]),"EP")

        if i==1: graph.Draw("AELPZ")
        elif i>1: graph.Draw("ELPZ same")

    leg.Draw()       
   
    c.SetLogy();
    for i,graph in enumerate(graphs_norm):
        graph.SetMinimum(0.03)
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6.6)
        if laser and recovery: graph.GetXaxis().SetLimits(0,40)
        #if "far" in name:
         #   graph.SetMinimum(0.15)
        graph.SetMaximum(3)
        #if "v9" not in outname:
         #   if i==0: graph.Draw("AELPZ")
          #  else: graph.Draw("ELPZ same")
        #else:
        if not ifb and ("mat" in legtype or variation >= 19):
            graph.SetLineStyle(2)
            if "Ref" in graphlabel[i]: continue 
            f1 = ROOT.TF1("f"+str(i),"expo",0,6.3)
            f1.SetLineColor(colors[i])
            f2 = ROOT.TF1("f2"+str(i),"expo",0,6.3)
            f2.SetLineColor(colors[i])
            f2.SetLineStyle(7)
            #f1.FixParameter(0,0.)
            #f1.FixParameter(4,0.)
            gclone = graph.Clone()
            gclone.Fit("f2"+str(i),"R")

            #if i==1: gclone.Draw("AELPZ")
            #elif i>1:  gclone.Draw("ELPZ same")
            
           # if i>0: gclone.Draw("ELPZ same")
           
            print "All points, fit parameter = "+str(f2.GetParameter(1))+" +/- "+str(f2.GetParError(1))
            doseconstant_allpoints=0
            if abs(f2.GetParameter(1)) > 0:
                doseconstant_allpoints= -1./f2.GetParameter(1)
                
            graph.RemovePoint(0)
            graph.Fit("f"+str(i),"R")
            f1.Draw("same")
            print "Remove first point, fit parameter = "+str(f1.GetParameter(1))+" +/- "+str(f1.GetParError(1))
            doseconstant=0
            doseconstant_unc=0
            if abs(f1.GetParameter(1)) > 0:
                doseconstant= -1./f1.GetParameter(1)
                doseconstant_unc= abs(doseconstant-doseconstant_allpoints)
                #doseconstant_unc= (0.2 + abs(f1.GetParError(1)/f1.GetParameter(1)) )/abs(f1.GetParameter(1))
                print "Dose constant = "+str(doseconstant)+" +/- "+str(doseconstant_unc)

            ##store dose constant in text file
            ##open old text file database. and copy through to new text file, looking for line corresponding to this channel
            ## if old line is found, replace it with new line
            ## else, add new line at end.
            
            outFile = open("DoseConstants_new.txt","w")
            try:
                with open("DoseConstants.txt") as channelrows:
                    for channelrow in channelrows:
                        if str(channels[i]) not in channelrow:
                            outFile.write(channelrow)
            except:
                print "Starting new DoseConstants.txt file"
            print "Adding channel ",str(channels[i])
            dose_rate = 15.9*0.25*dose_map[positions[i]]
            if "Ref" not in positions[i]: 
                beamdist = dist_map[positions[i]]
                outFile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s" %(graphlabel[i],str(channels[i]),beamdist,str(doseconstant),str(doseconstant_unc),str(dose_rate),str(0.15*dose_rate)+"\n"))                    
            
            try: os.remove("DoseConstants.txt")
            except: print "Creating new file"
            os.rename("DoseConstants_new.txt", "DoseConstants.txt")

        
        if i==1: graph.Draw("AELPZ")
        elif i>1: graph.Draw("ELPZ same")
                
    leg.Draw()
    if recovery:
        c.Print("doseplots/log_recovery_"+outname+"_"+name+"_vs_days_relative_change.pdf")
    else:
        if ifb: c.Print("doseplots/log_"+outname+"_"+name+"_vs_lumi_relative_change.pdf")
        else: c.Print("doseplots/log_"+outname+"_"+name+"_vs_dose_relative_change.pdf")
    c.Close()


    
      
print v_mean

     
