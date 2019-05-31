"""
*  author: carlos p
*  purpose: plot HMO r(t) vs t
*
"""

print("read_HMO_r_t.py")

''' controls '''
nlines=1
tId_1=1   # [0.1,1,10]

aId_0=0   # [0,0.1,1,10]
aId_1=1   # [0,0.1,1,10]
aId_2=2   # [0,0.1,1,10]

dtId_0=0
dtId_1=35
dtId_2=39

import sys
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from ROOT import TH1F,TStyle,kYellow

file=TFile("{}".format(sys.argv[1]), "read")
#file=TFile("data/hmo/HMO_dat_r1_v0_nstep1000_ndt40_nalpha4.root", "read")
#file=TFile("data/hmo/HMO_dat_r0_v1_nstep10000_ndt40_nalpha4.root", "read")
tree=file.Get("timeDat")
tree2=file.Get("freqDat")
tree3=file.Get("energyDat")

#tree.Scan("r:t:dt","","")
tree.GetEntry()
tree2.GetEntry()
tree3.GetEntry()

'''=====================================
    canvas 1
    ===================================='''
c1=TCanvas("c2", "canvas_2", 850, 400)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.1)
c1.SetRightMargin(0.1)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.1)

'''=====================================
    canvas 2
    ===================================='''
c2=TCanvas("c1", "canvas_1", 850, 400)
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderSize(0);
c2.SetLeftMargin(0.1)
c2.SetRightMargin(0.1)
c2.SetBottomMargin(0.15)
c2.SetTopMargin(0.1)

'''=====================================
    pads
    ===================================='''
pad1=TPad("pad1","xi_1_histogram",0,0.0,0.33,0.98);
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.1)
pad1.SetBottomMargin(0.15)
pad1.SetTopMargin(0.1)
pad1.Draw()
legend1 = TLegend(0.6, 0.70, .95, 0.87)
legend1.SetTextSize(0.05)

pad2=TPad("pad2","xi_2_histogram",0.33,0.0,0.66,0.98)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.1)
pad2.SetBottomMargin(0.15)
pad2.SetTopMargin(0.1)
pad2.Draw()
legend2 = TLegend(0.6, 0.70, .95, 0.87)

pad3=TPad("pad3","xi_3_histogram",0.66,0.0,0.99,0.98)
pad3.SetLeftMargin(0.15)
pad3.SetRightMargin(0.1)
pad3.SetBottomMargin(0.15)
pad3.SetTopMargin(0.1)
pad3.Draw()
legend3 = TLegend(0.6, 0.70, .95, 0.87)

'''=====================================
    canvas 1: plots
    ===================================='''
c1.cd()
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_0),"") # for pad 1
w0=tree.GetV3()[0]
a0=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_0),"same")
w1=tree.GetV3()[0]
a1=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_0),"same")
w2=tree.GetV3()[0]
a2=tree.GetV4()[0]

tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_1),"same") # for pad 2
ww0=tree.GetV3()[0]
aa0=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_1),"same")
ww1=tree.GetV3()[0]
aa1=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_1),"same")
ww2=tree.GetV3()[0]
aa2=tree.GetV4()[0]

tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_2),"same") # for pad 3
www0=tree.GetV3()[0]
aaa0=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_2),"same")
www1=tree.GetV3()[0]
aaa1=tree.GetV4()[0]
tree.Draw("r:t:Wdt:alpha","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_2),"same")
www2=tree.GetV3()[0]
aaa2=tree.GetV4()[0]

c1.Modified()
c1.Update()


'''=====================================
    canvas 2: pad 1 plot r and t at
                varying dt, alpha=0
    ===================================='''
c2.cd()
pad1.cd()

if nlines>=1:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_0),"") #alpha=0.1
    g0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g0.SetName("g0")
    g0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,a0))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(2)
    g0.GetXaxis().SetTitle('time')
    g0.GetYaxis().SetTitle('r(t)')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")
    legend1.AddEntry("g0","#Omega_0 dt=%.2f" % (w0),"l")

if nlines>=2:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_0),"same")
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(2)
    g1.Draw("L PLC PFC")
    legend1.AddEntry("g1","#Omega_{0} dt=%.2f" % (w1),"l")

if nlines>=3:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_0),"same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g2.SetName("g2")
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(2)
    g2.Draw("L PLC PFC")

#legend1.AddEntry("g2","#Omega_0 dt=%.2f" % (w2),"l")
legend1.Draw()


'''=====================================
    canvas 2: pad 2 plot r and t at
    varying dt, alpha=1
    ===================================='''
pad2.cd()
legend2.SetTextSize(0.05)

if nlines>=1:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_1),"") #alpha=0.1
    gg0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    gg0.SetName("gg0")
    gg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,aa0))
    gg0.SetMarkerStyle(2)
    gg0.SetLineWidth(2)
    gg0.GetXaxis().SetTitle('time')
    gg0.GetYaxis().SetTitle('r(t)')
    gg0.GetXaxis().SetTitleSize(0.06)
    gg0.GetYaxis().SetTitleSize(0.06)
    gg0.Draw("LA PLC PFC")
    legend2.AddEntry("gg0","#Omega_0 dt=%.2f" % (ww0),"l")

if nlines>=2:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_1),"same")
    gg1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    gg1.SetName("gg1")
    gg1.SetMarkerStyle(2)
    gg1.SetLineWidth(2)
    gg1.Draw("L PLC PFC")
    legend2.AddEntry("gg1","#Omega_0 dt=%.2f" % (ww1),"l")

if nlines>=3:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_1),"same")
    gg2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    gg2.SetName("gg2")
    gg2.SetMarkerStyle(2)
    gg2.SetLineWidth(2)
    gg2.Draw("L PLC PFC")
    legend2.AddEntry("gg2","#Omega_0 dt=%.2f" % (ww2),"l")

legend2.Draw()


'''=====================================
    canvas 2: pad 3 plot r and t at
    varying dt, alpha=1=2
    ===================================='''
pad3.cd()
legend3.SetTextSize(0.05)

if nlines>=1:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_0,aId_2),"") #alpha=0.1
    ggg0=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    ggg0.SetName("ggg0")
    ggg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.0f" % (tree.r0,tree.v0,aaa0))
    ggg0.SetMarkerStyle(2)
    ggg0.SetLineWidth(2)
    ggg0.GetXaxis().SetTitle('time')
    ggg0.GetYaxis().SetTitle('r(t)')
    ggg0.GetXaxis().SetTitleSize(0.06)
    ggg0.GetYaxis().SetTitleSize(0.06)
    ggg0.Draw("LA PLC PFC")
    legend3.AddEntry("ggg0","#Omega_0 dt=%.2f" % (www0),"l")

if nlines>=2:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_1,aId_2),"same")
    ggg1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    ggg1.SetName("ggg1")
    ggg1.SetMarkerStyle(2)
    ggg1.SetLineWidth(2)
    ggg1.Draw("L PLC PFC")
    legend3.AddEntry("ggg1","#Omega_0 dt=%.2f" % (www1),"l")

if nlines>=3:
    tree.Draw("r:t","temp_id=={} && dt_id=={} && alpha_id=={}".format(tId_1,dtId_2,aId_2),"same")
    ggg2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    ggg2.SetName("ggg2")
    ggg2.SetMarkerStyle(2)
    ggg2.SetLineWidth(2)
    ggg2.Draw("L PLC PFC")
    legend3.AddEntry("ggg2","#Omega_0 dt=%.2f" % (www2),"l")

legend3.Draw()

'''=====================================
    Settings
    ===================================='''
gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

pad1.SetGrid()
pad2.SetGrid()
pad3.SetGrid()

#legend.Draw()
c1.Modified()
c1.Update()

text=raw_input()


c2.Print("images/hmo/HMO_r{:.0f}_t_v{:.0f}_n10000000.pdf".format(tree.r0,tree.v0))
