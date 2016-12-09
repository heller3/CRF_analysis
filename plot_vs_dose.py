#!/usr/bin/env python
import os, sys, re
import ROOT
import os
from array import array

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.gStyle.SetLabelFont(42,"xyz")
ROOT.gStyle.SetLabelSize(0.04,"xyz")
#ROOT.gStyle.SetTitleFont(42)
ROOT.gStyle.SetTitleFont(42,"xyz")
ROOT.gStyle.SetTitleFont(42,"t")
#ROOT.gStyle.SetTitleSize(0.05)
ROOT.gStyle.SetTitleSize(0.05,"xyz")
ROOT.gStyle.SetPadBottomMargin(0.14)
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
shuntval = 3.87

#to apply a cut to runs considered (not ironed out yet)
applycut=False


if len(sys.argv) < 4:
    sys.exit('Please provide selection number, laser or ped mode, and xaxis mode') 

else:
    variation = int(sys.argv[1])
    if "ped" in sys.argv[2]: laser = False
    if "Mrad" in sys.argv[3]: ifb = False

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
    histname = "h1_energy_ieta"
    runlist = [280662,281740,282046,282368,282369,282779,282933,283403,283512,283660,283803, 283842, 283918, 283973, 284127]
    lumi =    [0.    ,1.1   ,2.31   ,2.70  ,2.72,  3.29, 4.01  ,5.90  ,6.8  ,7.0   ,7.24  , 7.66  , 8.81  , 9.53  , 10.0 ]
    #lumi =    [0.    ,1.1   ,2.1   ,2.49  ,2.51  ,  3.  ,3.73  ,5.61  ,6.52  ,6.73  ,6.96  , 7.37  , 8.53  , 9.25] original lumi calculation

#pedestals
else:
    histname = "h1_energy_const_ieta"
    runlist = [280646,281734,282040,282314,282772,282927,283321,283397,283522, 283662,283807, 283846, 283977, 284123]
    lumi =    [0.    ,1.1   ,2.31   ,2.71,  3.29, 4.01  , 5.44 ,5.90  ,6.8  ,7.0   ,7.24  , 7.66  , 9.53  , 10.0 ]
    #lumi =    [0.    ,1.1   ,2.1   ,2.5 ,  3.,    3.73  ,5.16  ,5.61  , 6.52  , 6.73  ,6.96  , 7.37  , 9.25]
    shunted = [0     ,0     ,1     ,1    ,1  ,     1    ,1     ,1     , 1      ,1      ,0     , 0    , 0] # indicates whether shunt was activated for run or not


legtype = "dist"    
outname = histname.strip('_ieta').replace("h1_","")
if fixshunt: outname = outname+"_shunt_corr"
if applycut and laser: outname = outname+"_cut"
        
if len(runlist) > 10:
    outname = outname+"_v11_"


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
        uf = float(hist.GetBinContent(0))/hist.GetEntries()
        if "const" in histname: ### this indicates pedestal run, divide by 8 to convert to single TS value
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
    else: max = 300
    graph.SetMaximum(max)
    graph.SetLineColor(colors[i])
    graph.SetMarkerColor(colors[i])
    graph.SetMarkerSize(1)
    graph.SetMarkerStyle(20)
    if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Average charge [fC]")
    else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Average charge [fC]")
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
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6)
        graph.SetMinimum(0.02)
        graph.SetMaximum(20)
        graph.SetLineColor(colors[i])
        graph.SetMarkerColor(colors[i])
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Average charge / reference charge")
        else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Average charge / reference charge")     #   graph.SetTitleFont(42);
        #leg.AddEntry(graph, "Chan. "+str(channels[i])+", position: "+positions[i],"EP")
        if "dist" in legtype:
            if "Ref" not in positions[i]: leg.AddEntry(graph, "R = "+dist_map[positions[i]]+"  "+str(channels[i]),"EP")
        else:
            if i>0: leg.AddEntry(graph, graphlabel[i]+" "+str(channels[i]),"EP")
                
        if i==1: graph.Draw("AELPZ")
        elif i>1: graph.Draw("ELPZ same")

    leg.Draw()       
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
        leg = ROOT.TLegend(0.47,leg_y_start,0.89,0.89)
    else:
        leg = ROOT.TLegend(0.65,0.65,0.89,0.89)
    leg.SetMargin(0.15)
    for i,graph in enumerate(graphs_norm):
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6.)
        graph.SetMinimum(0)
        graph.SetMaximum(1.15*max)
        graph.SetLineColor(colors[i])
        graph.SetMarkerColor(colors[i])
        graph.SetMarkerSize(1)
        graph.SetMarkerStyle(20)
        if ifb: graph.SetTitle(title+"; Integrated Luminosity [fb^{-1}]; Remaining relative brightness")
        else: graph.SetTitle(title+"; Integrated Dose [Mrad]; Remaining relative brightness")
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
        graph.SetMinimum(0.07)
        if not ifb: graph.GetXaxis().SetLimits(-0.6,6.)
        #if "far" in name:
         #   graph.SetMinimum(0.15)
        graph.SetMaximum(3)
        #if "v9" not in outname:
         #   if i==0: graph.Draw("AELPZ")
          #  else: graph.Draw("ELPZ same")
        #else:
        if i==1: graph.Draw("AELPZ")
        elif i>1: graph.Draw("ELPZ same")
                
    leg.Draw()
    if ifb: c.Print("doseplots/log_"+outname+"_"+name+"_vs_lumi_relative_change.pdf")
    else: c.Print("doseplots/log_"+outname+"_"+name+"_vs_dose_relative_change.pdf")
    c.Close()


    
      
print v_mean

     