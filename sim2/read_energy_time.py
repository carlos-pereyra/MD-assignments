"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import numpy as np
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall

nlines=2

file=TFile("data/particle_in_box_data.root", "read")
tree=file.Get("time_sequence")
tree2=file.Get("system_energy")
tree3=file.Get("avg_system_energy")
tree4=file.Get("system")

tree.GetEntry()
tree2.GetEntry()
tree3.GetEntry()
tree4.GetEntry()

#tree.Scan("x:y:ke:u_lj:u_wall","id==9")    #local issues?
tree2.Scan("u_wall_avg:u_lj_avg:dt:t","")  #energy vs time
#tree3.Scan("","")                          #energy vs dt
#tree4.Scan()

''' Canvas '''
c1=TCanvas("c1", "canvas 1", 600, 500)
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderMode(0);
c1.GetFrame().SetBorderSize(0);
c1.SetLeftMargin(0.15)
c1.SetRightMargin(0.15)
c1.SetBottomMargin(0.15)
c1.SetTopMargin(0.15)
legend = TLegend(0.50, 0.83, .95, 0.99)

if nlines>=1:
    tree2.Draw("ke_avg:t","t<20 && 0.009<dt<0.011","")
    g1=TGraph(tree2.GetSelectedRows(),tree2.GetV2(),tree2.GetV1())
    g1.SetName("g1")
    '''g1.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
                                                            ,tree.box.atoms[0].ry0
                                                            ,tree.box.atoms[0].vx0
                                                            ,tree.box.atoms[0].vy0))'''
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(3)
    g1.GetXaxis().SetTitle('Time')
    g1.GetYaxis().SetTitle('wall force_x')
    g1.GetXaxis().SetTitleSize(0.06)
    g1.GetYaxis().SetTitleSize(0.06)
    g1.Draw("LA PLC PFC")

if nlines>=2:
    tree2.Draw("ke_avg:t","t<20 && 0.019<dt<0.021","same")
    g2=TGraph(tree2.GetSelectedRows(),tree2.GetV2(),tree2.GetV1())
    g2.SetName("g2")
    '''g1.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
        ,tree.box.atoms[0].ry0
        ,tree.box.atoms[0].vx0
        ,tree.box.atoms[0].vy0))'''
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(3)
    g2.Draw("L PLC PFC")

#gStyle.SetOptStat(0)
gStyle.SetPalette(kRainBow)
#gStyle.SetOptTitle(0)
#gROOT.SetStyle("ATLAS")

c1.Modified()
c1.Update()

'''
gPad.Print("images/fx_x_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0
                                                                    ,tree.box.atoms[0].ry0
                                                                    ,tree.box.atoms[0].vx0
                                                                    ,tree.box.atoms[0].vy0))
'''

'''
#canvas 2
'''
'''
#c2=TCanvas("c2", "canvas 2", 600, 500)
#c2.GetFrame().SetBorderMode(0);
#c2.GetFrame().SetBorderMode(0);
#c2.GetFrame().SetBorderSize(0);
#c2.SetLeftMargin(0.15)
#c2.SetRightMargin(0.15)
#c2.SetBottomMargin(0.15)
#c2.SetTopMargin(0.06)
#c2.cd()

#graph2
tree.Draw("vx:t","id==0 && t<20","same")
g2=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g2.SetName("g2")
g2.SetLineWidth(4)
g2.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
                                                        ,tree.box.atoms[0].ry0
                                                        ,tree.box.atoms[0].vx0
                                                        ,tree.box.atoms[0].vy0))
g2.Draw("L PLC PFC")


#graph3
tree.Draw("box.atoms.r.y:box.atoms.t","id==0 && t<20","same")
g3=TGraph(tree.GetSelectedRows(),tree.GetV2(), tree.GetV1())
g3.SetName("g3")
g3.SetLineWidth(4)
g3.GetXaxis().SetTitleSize(0.06)
g3.GetYaxis().SetTitleSize(0.06)
g3.SetTitle('(x0,y0)=({},{}), (vx0,vy0)=({},{})'.format(tree.box.atoms[0].rx0
                                                        ,tree.box.atoms[0].ry0
                                                        ,tree.box.atoms[0].vx0
                                                        ,tree.box.atoms[0].vy0))
g3.Draw("L PLC PFC")

#output
#gROOT.SetStyle("ATLAS");
legend.SetTextSize(0.06)
legend.AddEntry("g1","U_wall vs time","l")
legend.AddEntry("g2","#xi vs time","l")
legend.AddEntry("g3","#nu vs time","l")
legend.Draw()
c.Modified()
c.Update()


gPad.Print("images/energy_wall_x{:.0f}_y{:.0f}_vx{:.0f}_vy{:.0f}.pdf".format(tree.box.atoms[0].rx0
                                                                        ,tree.box.atoms[0].ry0
                                                                        ,tree.box.atoms[0].vx0
                                                                        ,tree.box.atoms[0].vy0))'''


text=raw_input()
