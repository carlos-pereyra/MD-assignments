""" ============================================================
*  author: carlos p
*  purpose: plot average energy as a function of time step (dt)
*
============================================================ """

print("read_e_dt.py")
atom=0

import sys
import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall

file=TFile("{}".format(sys.argv[1]), "read")
#file=TFile("data/AtomBox_1atoms_1000tSteps_20dts_vx0_2_vy0_3.root", "read")
tree=file.Get("timeSequence")
tree2=file.Get("averagedEnergy")
tree3=file.Get("avgSystemEnergy")
tree4=file.Get("systemParameters")

tree.GetEntry()
tree2.GetEntry()
tree3.GetEntry()
tree4.GetEntry()

#tree.Scan("x:y:ke:u_lj:u_wall","id==9")    #local issues?
#tree2.Scan("u_wall_avg:u_lj_avg:dt:t","")  #energy vs time
#tree3.Scan("","")                          #energy vs dt
#tree4.Scan()

'''======================== Canvas 1 ==========================='''
c1=TCanvas("c1", "canvas 1", 600, 500)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.15)
legend = TLegend(0.50, 0.83, .95, 0.99)
legend.SetTextSize(0.04)

#graph 1
tree3.Draw("total_ke_v_avg:dt","total_ke_v_avg<10","")
g1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
g1.SetName("g1")
g1.SetTitle("vx_{0}=%d vy_{0}=%d" % (tree4.vx0,tree4.vy0))
g1.SetMarkerStyle(2)
g1.SetLineWidth(3)
g1.GetXaxis().SetTitle('Time')
#g1.GetYaxis().SetTitle('wall force_x')
g1.GetXaxis().SetTitleSize(0.06)
g1.GetYaxis().SetTitleSize(0.06)
g1.Draw("LA PLC PFC")
legend.AddEntry("g1","average ke w/ velocity=v","l")

tree3.Draw("total_e_v_avg:dt","","same")
gg1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
gg1.SetName("gg1")
gg1.SetMarkerStyle(2)
gg1.SetLineWidth(3)
gg1.Draw("L PLC PFC")
legend.AddEntry("gg1","average e w/ velocity=v","l")

'''=====================================
    Settings
    ===================================='''
gStyle.SetPalette(kRainBow)
#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0)
#gROOT.SetStyle("ATLAS")
legend.Draw()

c1.SetGrid()
c1.Modified()
c1.Update()

text=raw_input()

'''
gPad.Print("images/energyVtime_{}nAtoms_{}tSteps_{}dtSteps.pdf".format(tree4.natom,tree4.ntime,tree4.ndt))
'''
