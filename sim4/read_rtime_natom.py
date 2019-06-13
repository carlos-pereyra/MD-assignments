""" ============================================================
*  author: carlos p
*  purpose: plot potential and force
*
============================================================ """

print("read_p_natom.py")
DBG=0
nlines=5
dt_id=0
ylimit=10
import sys
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend,kAurora,kRainBow,kRust,kFall
from ROOT import TPaveText,TStyle,TMultiGraph,TGraphErrors

file=TFile("{}".format(sys.argv[1]), "read")
tree=file.Get("averageEnergy")
#tree.Scan("Sum$(pressure)/Length$(pressure):Sum$(natom)/Length$(natom)","","colsize=25")
if DBG: tree.Scan("runtime:natom","","colsize=25")
tree.GetEntry()

'''========================
    Canvas 1
    ==========================='''
c1=TCanvas("c", "canvas 1", 600, 600)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.10)
#legend = TLegend(0.60, 0.60, .98, 0.78)
#legend.SetTextSize(0.04)
c1.cd()

pad1=TPad("pad1","potential energy",0,0.06,0.99,0.99);
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.1)
pad1.SetBottomMargin(0.15)
pad1.SetTopMargin(0.1)
pad1.Draw()
legend1 = TLegend(0.5, 0.70, .95, 0.87)

pad2=TPad("pad2","force",0.0,0.0,0.99,0.05)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.1)
pad2.SetBottomMargin(0.15)
pad2.SetTopMargin(0.1)
pad2.Draw()
legend2 = TLegend(0.5, 0.40, .95, 0.57)

#mg=TMultiGraph()

'''========================
    Pad 1
    ==========================='''
pad1.cd()

tree.SetLineColor(1)
tree.SetTitle('normalized values')
tree.Draw("runtime:natom","","")
g0=TGraphErrors(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g0.SetName("g0")
#g0.SetTitle("Potential-Spline Diagram")
g0.SetMarkerStyle(2)
g0.SetMarkerColor(2)
g0.SetLineWidth(2)
g0.SetLineColor(1)
g0.SetLineStyle(1);
g0.GetXaxis().SetTitle('Natoms')
g0.GetYaxis().SetTitle('Elapsed Runtime')
g0.GetXaxis().SetTitleSize(0.06)
g0.GetYaxis().SetTitleSize(0.06)
#g0.GetXaxis().SetLimits(0.6,2)
#g0.GetYaxis().SetLimits(-10,10)
g0.Draw("ALP")
#legend1.AddEntry("g0","Try","lp")
#legend1.Draw()

#pt.AddText("#frac{2s}{#pi#alpha^{2}} ")
#pt.AddText("#vx_{0}=%d, vy_{0}=%d ")

'''
tree.Draw("ekuAvg:dt","","same")
g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g1.GetXaxis().SetLimits(0.6,2)
g1.SetMarkerStyle(9)
g1.SetLineWidth(2)
g1.SetLineColor(2)
g1.SetLineStyle(2);
g1.SetName("g1")
g1.Draw("L")
legend1.AddEntry("g1","<ek_{u}>","l")

tree.Draw("epvAvg:dt","","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g2.GetXaxis().SetLimits(0.6,2)
#g2.SetMarkerStyle(9)
g2.SetLineWidth(2)
g2.SetLineColor(1)
g2.SetLineStyle(1);
g2.SetName("g2")
g2.Draw("L")
legend1.AddEntry("g2","<ep_{v}>","l")

tree.Draw("epuAvg:dt","","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g3.GetXaxis().SetLimits(0.6,2)
#g3.SetMarkerStyle(2)
g3.SetLineWidth(2)
g3.SetLineColor(1)
g3.SetLineStyle(1);
g3.SetName("g3")
g3.Draw("L")
legend1.AddEntry("g3","<ep_{u}>","l")

legend1.Draw()'''

'''========================
    Pad 2
    ==========================='''

'''
pad2.cd()

#tree.SetLineColor(2)
tree.Draw("f_old:r","f<10","")
g4=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g4.GetXaxis().SetLimits(0.6,2)
#g4.SetMarkerStyle(2)
g4.SetLineWidth(2)
g4.SetLineColor(1)
g4.SetLineStyle(1);
g4.GetXaxis().SetTitle('x-direction')
g4.GetYaxis().SetTitle('normalized force')
g4.GetXaxis().SetTitleSize(0.06)
g4.GetYaxis().SetTitleSize(0.06)
g4.SetName("g4")
g4.Draw("LA")
legend2.AddEntry("g4","traditional lj force ","l")

tree.Draw("f:r","","same")
g5=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g5.GetXaxis().SetLimits(0.6,2)
#g5.SetMarkerStyle(2)
g5.SetLineWidth(2)
g5.SetLineColor(2)
g5.SetLineStyle(2);
g5.SetName("g5")
#g5.Draw("L")
#legend2.AddEntry("g5","spline lj force","l")

tree.Draw("fx:r","","same")
g6=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g6.GetXaxis().SetLimits(0.6,2)
#g6.SetMarkerStyle(2)
g6.SetLineWidth(2)
g6.SetLineColor(2)
g6.SetLineStyle(2);
g6.SetName("g6")
g6.Draw("L")
legend2.AddEntry("g6","spline lj force-x","l")

legend2.Draw()
'''

'''========================
    Settings
    ==========================='''
gStyle.SetOptStat(0)
#gStyle.SetPalette(kRainBow)
#pt.Draw()
gStyle.SetOptStat(0)
gStyle.SetOptTitle(0)

pad1.SetGrid()
pad2.SetGrid()
#mg.Add(g0)
#mg.Add(g1)
#mg.Draw("LA")


#c1.SetFillColor( 42 )
#c1.GetFrame().SetFillColor(21)
#c1.GetFrame().SetBorderSize(12)
#c1.SetGrid()
c1.Modified()
c1.Update()

#gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

text=raw_input()

'''
gPad.Print("images/Langevin_OpenBounds_particle_{}ntraces_{}nAtoms_{}tSteps_{}dtSteps.pdf"
           .format(nlines,tree4.natom,tree4.ntime,tree4.ndt))
'''
