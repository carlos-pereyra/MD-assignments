"""
*  author: carlos p
*  purpose: plot HMO eneergy vs Wdt
*
"""

nlines=2
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from ROOT import TH1F,TStyle,kYellow

file=TFile("data/hmo/HMO_wNoise_dat_r0_v1_n10000.root", "read")
tree=file.Get("time")
tree_freq=file.Get("freq")
tree_energy=file.Get("energy")

#tree.Scan("r:t:dt","","")
tree.GetEntry()
tree_freq.GetEntry()
tree_energy.GetEntry()

'''===================================== canvas 2 ===================================='''
c2=TCanvas("c2", "canvas_2", 850, 800)
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderSize(0);
c2.SetLeftMargin(0.1)
c2.SetRightMargin(0.1)
c2.SetBottomMargin(0.15)
c2.SetTopMargin(0.1)

c2.cd()
tree_energy.Draw("alpha","alpha_id==0","") # for pad 1
alpha0=tree_energy.GetV1()[0]

tree_energy.Draw("alpha","alpha_id==1","same") # for pad 2
aalpha0=tree_energy.GetV1()[0]

tree_energy.Draw("alpha","alpha_id==2","same") # for pad 3
aaalpha0=tree_energy.GetV1()[0]

tree_energy.Draw("alpha","alpha_id==3","same") # for pad 3
aaaalpha0=tree_energy.GetV1()[0]

c2.Modified()
c2.Update()

'''===================================== canvas 1 ===================================='''
c1=TCanvas("c1", "canvas_1", 850, 800)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.01)
c1.SetRightMargin(0.01)
c1.SetBottomMargin(0.01)
c1.SetTopMargin(0.01)

pad1=TPad("pad1","energy_vs_dt_for_a_specific_alpha",0,0.5,0.5,0.99);
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.1)
pad1.SetBottomMargin(0.15)
pad1.SetTopMargin(0.1)
pad1.Draw()
legend1 = TLegend(0.1, 0.70, .95, 0.87)

pad2=TPad("pad2","energy_vs_dt_for_a_specific_alpha",0.5,0.5,0.99,0.99)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.1)
pad2.SetBottomMargin(0.15)
pad2.SetTopMargin(0.1)
pad2.Draw()
legend2 = TLegend(0.1, 0.70, .95, 0.87)

pad3=TPad("pad3","energy_vs_dt_for_a_specific_alpha",0,0,0.50,0.50)
pad3.SetLeftMargin(0.15)
pad3.SetRightMargin(0.1)
pad3.SetBottomMargin(0.15)
pad3.SetTopMargin(0.1)
pad3.Draw()
legend3 = TLegend(0.1, 0.70, .95, 0.87)

pad4=TPad("pad4","energy_vs_dt_for_a_specific_alpha",0.5,0,0.99,0.50)
pad4.SetLeftMargin(0.15)
pad4.SetRightMargin(0.1)
pad4.SetBottomMargin(0.15)
pad4.SetTopMargin(0.1)
pad4.Draw()
legend4 = TLegend(0.1, 0.70, .95, 0.87)


'''===================================== pad 1 ===================================='''
pad1.cd()
legend1 = TLegend(0.15, 0.70, .95, 0.87)
if nlines>=1:
    tree_energy.Draw("totalE_avg_v:Wdt","alpha_id==0 && totalE_avg_v<5","")
    g0=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    g0.SetName("g0")
    g0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,alpha0))
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(2)
    g0.GetXaxis().SetTitle('#Omega_0 dt')
    g0.GetYaxis().SetTitle('Total Energy')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.Draw("LA PLC PFC")

if nlines>=2:
    tree_energy.Draw("totalE_avg_u:Wdt","alpha_id==0","same")
    g1=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(2)
    g1.Draw("L PLC PFC")

legend1.SetTextSize(0.05)
legend1.AddEntry("g0","full step velocity (avg energy)","l")
legend1.AddEntry("g1","half step velocity (avg energy)","l")
legend1.Draw()

'''===================================== pad 2 ===================================='''
pad2.cd()
legend2 = TLegend(0.15, 0.70, .95, 0.87)
if nlines>=1:
    tree_energy.Draw("totalE_avg_v:Wdt","alpha_id==1 && totalE_avg_v<5","")
    gg0=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    gg0.SetName("gg0")
    gg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,aalpha0))
    gg0.SetMarkerStyle(2)
    gg0.SetLineWidth(2)
    gg0.GetXaxis().SetTitle('#Omega_0 dt')
    gg0.GetYaxis().SetTitle('Total Energy')
    gg0.GetXaxis().SetTitleSize(0.06)
    gg0.GetYaxis().SetTitleSize(0.06)
    gg0.Draw("LA PLC PFC")

if nlines>=2:
    tree_energy.Draw("totalE_avg_u:Wdt","alpha_id==1","same")
    gg1=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    gg1.SetName("gg1")
    gg1.SetMarkerStyle(2)
    gg1.SetLineWidth(2)
    gg1.Draw("L PLC PFC")

legend2.SetTextSize(0.05)
legend2.AddEntry("gg0","full step velocity (avg energy)","l")
legend2.AddEntry("gg1","half step velocity (avg energy)","l")
legend2.Draw()

'''===================================== pad 3 ===================================='''
pad3.cd()
legend3 = TLegend(0.15, 0.70, .95, 0.87)
if nlines>=1:
    tree_energy.Draw("totalE_avg_v:Wdt","alpha_id==2 && totalE_avg_v<5","")
    ggg0=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    ggg0.SetName("ggg0")
    ggg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,aaalpha0))
    ggg0.SetMarkerStyle(2)
    ggg0.SetLineWidth(2)
    ggg0.GetXaxis().SetTitle('#Omega_0 dt')
    ggg0.GetYaxis().SetTitle('Total Energy')
    ggg0.GetXaxis().SetTitleSize(0.06)
    ggg0.GetYaxis().SetTitleSize(0.06)
    ggg0.Draw("LA PLC PFC")

if nlines>=2:
    tree_energy.Draw("totalE_avg_u:Wdt","alpha_id==2","same")
    ggg1=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    ggg1.SetName("ggg1")
    ggg1.SetMarkerStyle(2)
    ggg1.SetLineWidth(2)
    ggg1.Draw("L PLC PFC")

legend3.SetTextSize(0.05)
legend3.AddEntry("ggg0","full step velocity (avg energy)","l")
legend3.AddEntry("ggg1","half step velocity (avg energy)","l")
legend3.Draw()

'''===================================== pad 4 ===================================='''
pad4.cd()
legend4 = TLegend(0.15, 0.70, .95, 0.87)
if nlines>=1:
    tree_energy.Draw("totalE_avg_v:Wdt","alpha_id==3 && totalE_avg_v<5","")
    gggg0=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    gggg0.SetName("gggg0")
    gggg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree.r0,tree.v0,aaaalpha0))
    gggg0.SetMarkerStyle(2)
    gggg0.SetLineWidth(2)
    gggg0.GetXaxis().SetTitle('#Omega_0 dt')
    gggg0.GetYaxis().SetTitle('Total Energy')
    gggg0.GetXaxis().SetTitleSize(0.06)
    gggg0.GetYaxis().SetTitleSize(0.06)
    gggg0.Draw("LA PLC PFC")

if nlines>=2:
    tree_energy.Draw("totalE_avg_u:Wdt","alpha_id==3","same")
    gggg1=TGraph(tree_energy.GetSelectedRows(),tree_energy.GetV2(),tree_energy.GetV1())
    gggg1.SetName("gggg1")
    gggg1.SetMarkerStyle(2)
    gggg1.SetLineWidth(2)
    gggg1.Draw("L PLC PFC")

legend4.SetTextSize(0.05)
legend4.AddEntry("gggg0","full step velocity (avg energy)","l")
legend4.AddEntry("gggg1","half step velocity (avg energy)","l")
legend4.Draw()

'''================================================================================'''

gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0);
#gROOT.SetStyle("ATLAS");

pad1.SetGrid()
pad2.SetGrid()
pad3.SetGrid()
pad4.SetGrid()

c1.SetGrid()
c1.Modified()
c1.Update()

text=raw_input()

c1.Print("images/hmo_data/HMO_energy_wNoise_r%.0f_v%.0f_n%.0f.pdf" % (tree.r0,tree.v0,tree.size))
