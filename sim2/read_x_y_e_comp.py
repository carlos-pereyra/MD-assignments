""" ============================================================
*  author: carlos p
*  purpose: plot particle motion
*
============================================================ """

print("read_x_y_e_comp.py")

nlines=1
atom=0
import sys
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend,kAurora,kRainBow,kRust,kFall
from ROOT import TPaveText

file=TFile("{}".format(sys.argv[1]), "read")
#file=TFile("data/AtomBox_1atoms_1000tSteps_20dts.root", "read")
tree=file.Get("timeSequence")
tree2=file.Get("averagedEnergy")
tree3=file.Get("avgSystemEnergy")
tree4=file.Get("systemParameters")

#tree.Scan("x:y:ke:u_lj:u_wall","id==9")    #local issues?
#tree2.Scan("u_wall_avg:u_lj_avg:dt:t","")  #energy vs time
#tree3.Scan("","")                          #energy vs dt
#tree4.Scan()

tree.GetEntry()
tree2.GetEntry()
tree3.GetEntry()
tree4.GetEntry()

'''========================
    Canvas 1
    ==========================='''
c1=TCanvas("c", "canvas 1", 500, 900)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.10)

pad1=TPad("pad1","particle motion 1",0.00,0.66,0.98,0.99);
pad1.SetLeftMargin(0.15)
pad1.SetRightMargin(0.1)
pad1.SetBottomMargin(0.15)
pad1.SetTopMargin(0.1)
pad1.Draw()
legend1 = TLegend(0.60, 0.80, .98, 0.98)
legend1.SetTextSize(0.04)

pad2=TPad("pad2","particle energy 1",0.00,0.33,0.98,0.66)
pad2.SetLeftMargin(0.15)
pad2.SetRightMargin(0.1)
pad2.SetBottomMargin(0.15)
pad2.SetTopMargin(0.1)
pad2.Draw()
legend2 = TLegend(0.60, 0.80, .98, 0.98)

pad3=TPad("pad3","particle energy 1",0.00,0.00,0.98,0.33)
pad3.SetLeftMargin(0.15)
pad3.SetRightMargin(0.1)
pad3.SetBottomMargin(0.15)
pad3.SetTopMargin(0.1)
pad3.Draw()
legend3 = TLegend(0.60, 0.80, .98, 0.98)

'''=====================================
    pad 1
    ===================================='''
pad1.cd()

if nlines>=1 or atom==0:
    tree.Draw("y:x","Iteration$==0 && dt_id==0","")
    g1=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g1.SetName("g1")
    g1.SetTitle("Motion Trace, vx_{0}=%d vy_{0}=%d" % (tree4.vx0,tree4.vy0))
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(3)
    g1.GetXaxis().SetTitle('#xi')
    g1.GetYaxis().SetTitle('#nu')
    g1.GetXaxis().SetTitleSize(0.06)
    g1.GetYaxis().SetTitleSize(0.06)
    g1.GetXaxis().SetLimits(0,20)
    g1.GetYaxis().SetLimits(0,20)
    g1.Draw("LA PLC PFC")
    if atom==0: legend1.AddEntry("g1","Particle 0 dt={:.2e}".format(tree.dt),"l")

if nlines>=2 or atom==1:
    tree.Draw("y:x","Iteration$==1 && dt_id==0","same")
    g2=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g2.SetName("g2")
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(3)
    g2.Draw("L PLC PFC")
    if atom==1: legend1.AddEntry("g2","Particle 1 dt={:.2f}".format(tree.dt),"l")

if nlines>=3 or atom==2:
    tree.Draw("y:x","Iteration$==2 && dt_id==0","same")
    g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g3.SetName("g3")
    g3.SetMarkerStyle(2)
    g3.SetLineWidth(3)
    g3.Draw("L PLC PFC")
    if atom==2: legend1.AddEntry("g3","Particle 2 dt={:.2f}".format(tree.dt),"l")

if nlines>=4 or atom==3:
    tree.Draw("y:x","Iteration$==3 && dt_id==0","same")
    g4=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g4.SetName("g4")
    g4.SetMarkerStyle(2)
    g4.SetLineWidth(3)
    g4.Draw("L PLC PFC")
    if atom==3: legend1.AddEntry("g4","Particle 3 dt={:.2f}".format(tree.dt),"l")

if nlines>=5 or atom==4:
    tree.Draw("y:x","Iteration$==4 && dt_id==0","same")
    g5=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g5.SetName("g5")
    g5.SetMarkerStyle(2)
    g5.SetLineWidth(3)
    g5.Draw("L PLC PFC")
    if atom==4: legend1.AddEntry("g5","Particle 4 dt={:.2f}".format(tree.dt),"l")

if nlines>=6 or atom==5:
    tree.Draw("y:x","Iteration$==5 && dt_id==0","same")
    g6=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g6.SetName("g6")
    g6.SetMarkerStyle(2)
    g6.SetLineWidth(3)
    g6.Draw("L PLC PFC")
    if atom==5: legend1.AddEntry("g6","Particle 5 dt={:.2f}".format(tree.dt),"l")

if nlines>=7 or atom==6:
    tree.Draw("y:x","Iteration$==6 && dt_id==0","same")
    g7=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g7.SetName("g7")
    g7.SetMarkerStyle(2)
    g7.SetLineWidth(3)
    g7.Draw("L PLC PFC")
    if atom==6: legend1.AddEntry("g7","Particle 6 dt={:.2f}".format(tree.dt),"l")

if nlines>=8 or atom==7:
    tree.Draw("y:x","Iteration$==7 && dt_id==0","same")
    g8=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g8.SetName("g8")
    g8.SetMarkerStyle(2)
    g8.SetLineWidth(3)
    g8.Draw("L PLC PFC")
    if atom==7: legend1.AddEntry("g8","Particle 7 dt={:.2f}".format(tree.dt),"l")

if nlines>=9 or atom==8:
    tree.Draw("y:x","Iteration$==8 && dt_id==0","same")
    g9=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g9.SetName("g9")
    g9.SetMarkerStyle(2)
    g9.SetLineWidth(3)
    g9.Draw("L PLC PFC")
    if atom==8: legend1.AddEntry("g9","Particle 8 dt={:.2f}".format(tree.dt),"l")

if nlines>=10 or atom==9:
    tree.Draw("y:x","Iteration$==9 && dt_id==0","same")
    g10=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g10.SetName("g10")
    g10.SetMarkerStyle(2)
    g10.SetLineWidth(3)
    g10.Draw("L PLC PFC")
    if atom==9: legend1.AddEntry("g10","Particle 9 dt={:.2f}".format(tree.dt),"l")

if nlines>=11 or atom==10:
    tree.Draw("y:x","Iteration$==10 && dt_id==0","same")
    g11=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g11.SetName("g11")
    g11.SetMarkerStyle(2)
    g11.SetLineWidth(3)
    g11.Draw("L PLC PFC")
    if atom==10: legend1.AddEntry("g11","Particle 10 dt={:.2f}".format(tree.dt),"l")

if nlines>=12 or atom==11:
    tree.Draw("y:x","Iteration$==11 && dt_id==0","same")
    g12=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g12.SetName("g12")
    g12.SetMarkerStyle(2)
    g12.SetLineWidth(3)
    g12.Draw("L PLC PFC")
    if atom==11: legend1.AddEntry("g12","Particle 11 dt={:.2f}".format(tree.dt),"l")

if nlines>=13 or atom==12:
    tree.Draw("y:x","Iteration$==12 && dt_id==0","same")
    g13=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g13.SetName("g13")
    g13.SetMarkerStyle(2)
    g13.SetLineWidth(3)
    g13.Draw("L PLC PFC")
    if atom==12: legend1.AddEntry("g13","Particle 12 dt={:.2f}".format(tree.dt),"l")

if nlines>=14 or atom==13:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g14=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g14.SetName("g14")
    g14.SetMarkerStyle(2)
    g14.SetLineWidth(3)
    g14.Draw("L PLC PFC")
    if atom==13: legend1.AddEntry("g14","Particle 13 dt={:.2f}".format(tree.dt),"l")

if nlines>=15 or atom==14:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g15=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g15.SetName("g15")
    g15.SetMarkerStyle(2)
    g15.SetLineWidth(3)
    g15.Draw("L PLC PFC")
    if atom==14: legend1.AddEntry("g15","Particle 14 dt={:.2f}".format(tree.dt),"l")

if nlines>=16 or atom==15:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g16=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g16.SetName("g16")
    g16.SetMarkerStyle(2)
    g16.SetLineWidth(3)
    g16.Draw("L PLC PFC")
    if atom==15: legend1.AddEntry("g16","Particle 15 dt={:.2f}".format(tree.dt),"l")

if nlines>=17 or atom==16:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g17=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g17.SetName("g17")
    g17.SetMarkerStyle(2)
    g17.SetLineWidth(3)
    g17.Draw("L PLC PFC")
    if atom==16: legend1.AddEntry("g17","Particle 16 dt={:.2f}".format(tree.dt),"l")

if nlines>=18 or atom==17:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g18=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g18.SetName("g18")
    g18.SetMarkerStyle(2)
    g18.SetLineWidth(3)
    g18.Draw("L PLC PFC")
    if atom==17: legend1.AddEntry("g18","Particle 17 dt={:.2f}".format(tree.dt),"l")

if nlines>=19 or atom==18:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g19=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g19.SetName("g19")
    g19.SetMarkerStyle(2)
    g19.SetLineWidth(3)
    g19.Draw("L PLC PFC")
    if atom==18: legend1.AddEntry("g19","Particle 18 dt={:.2f}".format(tree.dt),"l")

if nlines>=20 or atom==19:
    tree.Draw("y:x","Iteration$==13 && dt_id==0","same")
    g20=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
    g20.SetName("g20")
    g20.SetMarkerStyle(2)
    g20.SetLineWidth(3)
    g20.Draw("L PLC PFC")
    if atom==19: legend1.AddEntry("g20","Particle 19 dt={:.2f}".format(tree.dt),"l")

#pt.AddText("#frac{2s}{#pi#alpha^{2}} ")
#pt.AddText("#vx_{0}=%d, vy_{0}=%d ")
legend1.Draw()

'''=====================================
    pad 2
    ===================================='''
pad2.cd()

tree.Draw("ke_u:x","Iteration$=={} && dt_id==0".format(atom),"")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g3.SetName("g3")
g3.SetTitle("Kinetic Energies, Particle={}".format(atom))
g3.SetMarkerStyle(2)
g3.SetLineWidth(3)
g3.GetXaxis().SetTitle('#xi')
g3.GetYaxis().SetTitle('KE')
g3.GetXaxis().SetTitleSize(0.06)
g3.GetYaxis().SetTitleSize(0.06)
g3.GetXaxis().SetLimits(0,20)
g3.GetYaxis().SetLimits(0,20)
g3.Draw("LA PLC PFC")
legend2.AddEntry("g3","Kinetic Energy V dt={:.2f}".format(tree.dt),"l")

tree.Draw("ke_v:x","Iteration$=={} && dt_id==0".format(atom),"same")
g4=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g4.SetName("g4")
g4.SetLineWidth(3)
g4.Draw("L PLC PFC")
#pt.AddText("#frac{2s}{#pi#alpha^{2}} ")
#pt.AddText("#vx_{0}=%d, vy_{0}=%d ")
legend2.SetTextSize(0.04)
legend2.AddEntry("g4","Kinetic Energy U dt={:.2f}".format(tree.dt),"l")
legend2.Draw()

'''=====================================
    pad 3
    ===================================='''
pad3.cd()
tree.Draw("u_wall:x","Iteration$=={} && dt_id==0".format(atom),"")
g5=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g5.SetName("g5")
g5.SetTitle("Potential Energies, Particle={}".format(atom))
g5.SetMarkerStyle(2)
g5.SetLineWidth(3)
g5.GetXaxis().SetTitle('#xi')
g5.GetYaxis().SetTitle('PE')
g5.GetXaxis().SetTitleSize(0.06)
g5.GetYaxis().SetTitleSize(0.06)
g5.GetXaxis().SetLimits(0,20)
g5.GetYaxis().SetLimits(0,20)
g5.Draw("LA PLC PFC")
legend3.AddEntry("g5","Wall Potential dt={:.2f}".format(tree.dt),"l")

tree.Draw("u_jones:x","Iteration$=={} && dt_id==0".format(atom),"same")
g6=TGraph(tree.GetSelectedRows(),tree.GetV2(),tree.GetV1())
g6.SetName("g6")
g6.SetLineWidth(3)
g6.Draw("L PLC PFC")
legend3.SetTextSize(0.04)
legend3.AddEntry("g6","Lennard Jones dt={:.2f}".format(tree.dt),"l")
legend3.Draw()

'''=====================================
    Settings
    ===================================='''
gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0)
#pt.Draw()

pad1.SetGrid()
pad2.SetGrid()
pad3.SetGrid()
c1.Modified()
c1.Update()

#gStyle.SetOptStat(0)
#gStyle.SetOptTitle(0);
gStyle.SetPalette(kRainBow);
#gROOT.SetStyle("ATLAS");

text=raw_input()

'''
gPad.Print("images/particle_{}ntraces_{}nAtoms_{}tSteps_{}dtSteps.pdf"
           .format(nlines,tree4.natom,tree4.ntime,tree4.ndt))
'''
