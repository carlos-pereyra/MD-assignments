"""
*  author: carlos p
*  purpose: plot energy vs dt
*
"""

print("read_e_dt_array.py")

nlines=1
import sys
import numpy as np
from array import array
import ROOT
from ROOT import TTree,TFile,gROOT,gStyle,gPad
from ROOT import TCanvas,TGraph,TPad,TBrowser
from ROOT import TLegend, kAurora, kRainBow, kRust, kFall
from ROOT import TH1F,TStyle,kYellow

file=TFile("{}".format(sys.argv[1]), "read")
if file.IsZombie():
    print("read_x_y.py error: cannont read file")
    exit(1)

#file=TFile("data/dat_5natom_1000nsteps_2vx_3vy_wrTime.root", "read")
tree=file.Get("timeSequence")
tree3=file.Get("avgSystemEnergy")
tree4=file.Get("systemParameters")
tree.GetEntry()
tree3.GetEntry()
tree4.GetEntry()
tree3.Scan("","temp_id==0 && alpha_id==0","")


'''=====================================
    canvas 2
    ===================================='''
c2=TCanvas("c2", "canvas_2", 850, 800)
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderMode(0);
c2.GetFrame().SetBorderSize(0);
c2.SetLeftMargin(0.1)
c2.SetRightMargin(0.1)
c2.SetBottomMargin(0.15)
c2.SetTopMargin(0.1)

c2.cd()
tree3.Draw("alpha:temp","alpha_id==0 && temp_id==0","") # for pad 1
alpha0=tree3.GetV1()[0]
temp0=tree3.GetV2()[0]

tree3.Draw("alpha:temp","alpha_id==1 && temp_id==1","same") # for pad 2
aalpha0=tree3.GetV1()[0]
ttemp0=tree3.GetV2()[0]

tree3.Draw("alpha:temp","alpha_id==2 && temp_id==2","same") # for pad 3
aaalpha0=tree3.GetV1()[0]
tttemp0=tree3.GetV2()[0]

tree3.Draw("alpha:temp","alpha_id==3 && temp_id==0","same") # for pad 3
aaaalpha0=tree3.GetV1()[0]
ttttemp0=tree3.GetV2()[0]

c2.Modified()
c2.Update()

'''=====================================
    canvas 1
    ===================================='''
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


'''=====================================
    pad 1: alpha 1
    ===================================='''
pad1.cd()
legend1 = TLegend(0.2, 0.75, .95, 0.85)
legend1.SetTextSize(0.05)

if nlines==1: #temp=1
    tree3.Draw("ke_v_navg:dt","temp_id==1 && alpha_id==0 && ke_v_navg<20","")
    g0=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g0.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,alpha0))
    g0.SetName("g0")
    g0.SetMarkerStyle(2)
    g0.SetLineWidth(2)
    g0.GetXaxis().SetTitle('dt')
    g0.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    g0.GetXaxis().SetTitleSize(0.06)
    g0.GetYaxis().SetTitleSize(0.06)
    g0.GetXaxis().SetLabelSize(0.05)
    g0.GetYaxis().SetLabelSize(0.05)
    g0.Draw("LA PLC PFC")
    legend1.AddEntry("g0","avg energy (vel=v) temp={}".format(ttemp0),"l")
    
    tree3.Draw("ke_u_navg:dt","temp_id==1 && alpha_id==0","same")
    g1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g1.SetName("g1")
    g1.SetMarkerStyle(2)
    g1.SetLineWidth(2)
    g1.Draw("L PLC PFC")
    legend1.AddEntry("g1","avg energy (vel=u) temp={:.1f}".format(ttemp0),"l")

if nlines==2: #temp=10
    tree3.Draw("ke_v_navg:dt","temp_id==2 && alpha_id==0 && ke_v_navg<30","same")
    g2=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g2.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,alpha0))
    g2.SetName("g2")
    g2.SetMarkerStyle(2)
    g2.SetLineWidth(2)
    g2.GetXaxis().SetTitle('dt')
    g2.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    g2.GetXaxis().SetTitleSize(0.06)
    g2.GetYaxis().SetTitleSize(0.06)
    g2.GetXaxis().SetLabelSize(0.05)
    g2.GetYaxis().SetLabelSize(0.05)
    #g2.GetXaxis().LabelsOption('h')
    g2.Draw("LA PLC PFC")
    legend1.AddEntry("g2","avg energy (vel=v) temp={:.1f}".format(tttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==2 && alpha_id==0","same")
    g3=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g3.SetName("g3")
    g3.SetMarkerStyle(2)
    g3.SetLineWidth(2)
    g3.Draw("L PLC PFC")
    legend1.AddEntry("g3","avg energy (vel=u) temp={:.1f}".format(tttemp0),"l")

if nlines==3: #temp=0.1
    tree3.Draw("ke_v_navg:dt","temp_id==0 && alpha_id==0 && ke_v_navg<15","same")
    g4=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g4.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,alpha0))
    g4.SetName("g4")
    g4.SetMarkerStyle(2)
    g4.SetLineWidth(2)
    g4.GetXaxis().SetTitle('dt')
    g4.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    g4.GetXaxis().SetTitleSize(0.06)
    g4.GetYaxis().SetTitleSize(0.06)
    g4.GetXaxis().SetLabelSize(0.05)
    g4.GetYaxis().SetLabelSize(0.05)
    g4.Draw("LA PLC PFC")
    legend1.AddEntry("g4","avg energy (vel=v) temp={:.1f}".format(temp0),"l")
    
    tree3.Draw("ke_u_navg:dt","temp_id==0 && alpha_id==0","same")
    g5=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    g5.SetName("g5")
    g5.SetMarkerStyle(2)
    g5.SetLineWidth(2)
    g5.Draw("L PLC PFC")
    legend1.AddEntry("g5","avg energy (vel=u) temp={:.1f}".format(temp0),"l")

legend1.Draw()

'''=====================================
    pad 2: alpha 0.1
    ===================================='''
pad2.cd()
legend2 = TLegend(0.2, 0.75, .95, 0.85)
legend2.SetTextSize(0.05)

if nlines==1: #temp=1
    tree3.Draw("ke_v_navg:dt","temp_id==1 && alpha_id==1 && ke_v_navg<20","")
    gg0=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg0.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aalpha0))
    gg0.SetName("gg0")
    gg0.SetMarkerStyle(2)
    gg0.SetLineWidth(2)
    gg0.GetXaxis().SetTitle('dt')
    gg0.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    gg0.GetXaxis().SetTitleSize(0.06)
    gg0.GetYaxis().SetTitleSize(0.06)
    gg0.GetXaxis().SetLabelSize(0.05)
    gg0.GetYaxis().SetLabelSize(0.05)
    gg0.Draw("LA PLC PFC")
    legend2.AddEntry("gg0","avg energy (vel=v) temp={:.1f}".format(ttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==1 && alpha_id==1","same")
    gg1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg1.SetName("gg1")
    gg1.SetMarkerStyle(2)
    gg1.SetLineWidth(2)
    gg1.Draw("L PLC PFC")
    legend2.AddEntry("gg1","avg energy (vel=u) temp={:.1f}".format(ttemp0),"l")

if nlines==2: #temp=10
    tree3.Draw("ke_v_navg:dt","temp_id==2 && alpha_id==1 && ke_v_navg<20","same")
    gg2=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg2.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aalpha0))
    gg2.SetName("gg2")
    gg2.SetMarkerStyle(2)
    gg2.SetLineWidth(2)
    gg2.GetXaxis().SetTitle('dt')
    gg2.GetYaxis().SetTitle('Total Average Energy (<e_{k}>)')
    gg2.GetXaxis().SetTitleSize(0.06)
    gg2.GetYaxis().SetTitleSize(0.06)
    gg2.GetXaxis().SetLabelSize(0.05)
    gg2.GetYaxis().SetLabelSize(0.05)
    gg2.Draw("LA PLC PFC")
    legend2.AddEntry("gg2","avg energy (vel=v) temp={:.1f}".format(tttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==2 && alpha_id==1","same")
    gg3=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg3.SetName("gg3")
    gg3.SetMarkerStyle(2)
    gg3.SetLineWidth(2)
    gg3.Draw("L PLC PFC")
    legend2.AddEntry("gg3","avg energy (vel=u) temp={:.1f}".format(tttemp0),"l")

if nlines==3: #temp=0.1
    tree3.Draw("ke_v_navg:dt","temp_id==0 && alpha_id==1 && ke_v_navg<15","same")
    gg4=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg4.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aalpha0))
    gg4.SetName("gg4")
    gg4.SetMarkerStyle(2)
    gg4.SetLineWidth(2)
    gg4.GetXaxis().SetTitle('dt')
    gg4.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    gg4.GetXaxis().SetTitleSize(0.06)
    gg4.GetYaxis().SetTitleSize(0.06)
    gg4.GetXaxis().SetLabelSize(0.05)
    gg4.GetYaxis().SetLabelSize(0.05)
    gg4.Draw("LA PLC PFC")
    legend2.AddEntry("gg2","avg energy (vel=v) temp={:.1f}".format(temp0),"l")
    
    tree3.Draw("ke_u_navg:dt","temp_id==0 && alpha_id==1","same")
    gg5=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gg5.SetName("gg5")
    gg5.SetMarkerStyle(2)
    gg5.SetLineWidth(2)
    gg5.Draw("L PLC PFC")
    legend2.AddEntry("gg2","avg energy (vel=u) temp={:.1f}".format(temp0),"l")

legend2.Draw()

'''=====================================
    pad 3: alpha 1
    ===================================='''
pad3.cd()
legend3 = TLegend(0.2, 0.75, .95, 0.85)
legend3.SetTextSize(0.05)

if nlines==1: #temp=1
    tree3.Draw("ke_v_navg:dt","temp_id==1 && alpha_id==2 && ke_v_navg<15","")
    ggg0=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg0.SetName("ggg0")
    ggg0.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aaalpha0))
    ggg0.SetMarkerStyle(2)
    ggg0.SetLineWidth(2)
    ggg0.GetXaxis().SetTitle('dt')
    ggg0.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    ggg0.GetXaxis().SetTitleSize(0.06)
    ggg0.GetYaxis().SetTitleSize(0.06)
    ggg0.GetXaxis().SetLabelSize(0.05)
    ggg0.GetYaxis().SetLabelSize(0.05)
    ggg0.Draw("LA PLC PFC")
    legend3.AddEntry("ggg0","avg energy (vel=v) temp={:.1f}".format(ttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==1 && alpha_id==2","same")
    ggg1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg1.SetName("ggg1")
    ggg1.SetMarkerStyle(2)
    ggg1.SetLineWidth(2)
    ggg1.Draw("L PLC PFC")
    legend3.AddEntry("ggg1","avg energy (vel=u) temp={:.1f}".format(ttemp0),"l")

if nlines==2: #temp=10
    tree3.Draw("ke_v_navg:dt","temp_id==2 && alpha_id==2 && ke_v_navg<15","same")
    ggg2=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg2.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aaalpha0))
    ggg2.SetName("ggg2")
    ggg2.SetMarkerStyle(2)
    ggg2.SetLineWidth(2)
    ggg2.GetXaxis().SetTitle('dt')
    ggg2.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    ggg2.GetXaxis().SetTitleSize(0.06)
    ggg2.GetYaxis().SetTitleSize(0.06)
    ggg2.GetXaxis().SetLabelSize(0.05)
    ggg2.GetYaxis().SetLabelSize(0.05)
    ggg2.Draw("LA PLC PFC")
    legend3.AddEntry("ggg2","avg energy (vel=v) temp={:.1f}".format(tttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==2 && alpha_id==2","same")
    ggg3=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg3.SetName("ggg3")
    ggg3.SetMarkerStyle(2)
    ggg3.SetLineWidth(2)
    ggg3.Draw("L PLC PFC")
    legend3.AddEntry("ggg3","avg energy (vel=u) temp={:.1f}".format(tttemp0),"l")

if nlines==3: #temp=0.1
    tree3.Draw("ke_v_navg:dt","temp_id==0 && alpha_id==2 && ke_v_navg<15","same")
    ggg4=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg4.SetTitle("vx_{0}=%.0f -- vy_{0}=%.0f -- #alpha=%.1f" % (tree4.vx0,tree4.vy0,aaalpha0))
    ggg4.SetName("ggg4")
    ggg4.SetMarkerStyle(2)
    ggg4.SetLineWidth(2)
    ggg4.GetXaxis().SetLabelSize(0.05)
    ggg4.GetYaxis().SetLabelSize(0.05)
    ggg4.GetXaxis().SetTitle('dt')
    ggg4.GetYaxis().SetTitle('Average Energy (<e_{k}>)')
    ggg4.GetXaxis().SetTitleSize(0.06)
    ggg4.GetYaxis().SetTitleSize(0.06)
    ggg4.Draw("LA PLC PFC")
    legend3.AddEntry("ggg4","avg energy (vel=v) temp={:.1f}".format(temp0),"l")
    
    tree3.Draw("ke_u_navg:dt","temp_id==0 && alpha_id==2","same")
    ggg5=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    ggg5.SetName("ggg5")
    ggg5.SetMarkerStyle(2)
    ggg5.SetLineWidth(2)
    ggg5.Draw("L PLC PFC")
    legend3.AddEntry("ggg5","avg energy (vel=u) temp={:.1f}".format(temp0),"l")

legend3.Draw()

'''=====================================
    pad 4: alpha 10
    ===================================='''
'''
pad4.cd()
legend4 = TLegend(0.2, 0.75, .95, 0.85)
legend4.SetTextSize(0.05)

if nlines==1: #temp 1
    tree3.Draw("ke_v_navg:dt","temp_id==1 && alpha_id==3 && ke_v_navg<20","")
    gggg0=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg0.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree4.vx0,tree4.vy0,aaaalpha0))
    gggg0.SetName("gggg0")
    gggg0.SetMarkerStyle(2)
    gggg0.SetLineWidth(2)
    gggg0.GetXaxis().SetTitle('dt')
    gggg0.GetYaxis().SetTitle('Total Average Energy')
    gggg0.GetXaxis().SetTitleSize(0.06)
    gggg0.GetYaxis().SetTitleSize(0.06)
    gggg0.GetXaxis().SetLabelSize(0.05)
    gggg0.GetYaxis().SetLabelSize(0.05)
    gggg0.Draw("LA PLC PFC")
    legend4.AddEntry("gggg0","avg energy (vel=v) temp={:.1f}".format(ttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==1 && alpha_id==3","same")
    gggg1=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg1.SetName("gggg1")
    gggg1.SetMarkerStyle(2)
    gggg1.SetLineWidth(2)
    gggg1.Draw("L PLC PFC")
    legend4.AddEntry("gggg1","avg energy (vel=u) temp={:.1f}".format(ttemp0),"l")

if nlines==2: #temp=10
    tree3.Draw("ke_v_navg:dt","temp_id==2 && alpha_id==3 && ke_v_navg<20","")
    gggg2=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg2.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree4.vx0,tree4.vx0,aaaalpha0))
    gggg2.SetName("gggg2")
    gggg2.SetMarkerStyle(2)
    gggg2.SetLineWidth(2)
    gggg2.GetXaxis().SetTitle('dt')
    gggg2.GetYaxis().SetTitle('Total Average Energy')
    gggg2.GetXaxis().SetTitleSize(0.06)
    gggg2.GetYaxis().SetTitleSize(0.06)
    gggg2.GetXaxis().SetLabelSize(0.05)
    gggg2.GetYaxis().SetLabelSize(0.05)
    gggg2.Draw("LA PLC PFC")
    legend4.AddEntry("gggg2","avg energy (vel=v) temp={:.1f}".format(tttemp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==2 && alpha_id==3","same")
    gggg3=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg3.SetName("gggg3")
    gggg3.SetMarkerStyle(2)
    gggg3.SetLineWidth(2)
    gggg3.Draw("L PLC PFC")
    legend4.AddEntry("gggg3","avg energy (vel=u) temp={:.1f}".format(tttemp0),"l")

if nlines==3: #temp=0.1
    tree3.Draw("ke_v_avg:dt","temp_id==0 && alpha_id==3 && ke_v_navg<20","")
    gggg4=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg4.SetTitle("r0=%.0f, v0=%.0f, alpha=%.1f" % (tree4.vx0,tree4.vy0,aaaalpha0))
    gggg4.SetName("gggg4")
    gggg4.SetMarkerStyle(2)
    gggg4.SetLineWidth(2)
    gggg4.GetXaxis().SetTitle('dt')
    gggg4.GetYaxis().SetTitle('Total Average Energy')
    gggg4.GetXaxis().SetTitleSize(0.06)
    gggg4.GetYaxis().SetTitleSize(0.06)
    gggg4.GetXaxis().SetLabelSize(0.05)
    gggg4.GetYaxis().SetLabelSize(0.05)
    gggg4.Draw("LA PLC PFC")
    legend4.AddEntry("gggg4","avg energy (vel=v) temp={:.1f}".format(temp0),"l")

    tree3.Draw("ke_u_navg:dt","temp_id==0 && alpha_id==3","same")
    gggg5=TGraph(tree3.GetSelectedRows(),tree3.GetV2(),tree3.GetV1())
    gggg5.SetName("gggg5")
    gggg5.SetMarkerStyle(2)
    gggg5.SetLineWidth(2)
    gggg5.Draw("L PLC PFC")
    legend4.AddEntry("gggg5","avg energy (vel=u) temp={:.1f}".format(temp0),"l")

legend4.Draw()
'''

'''=====================================
    Settings
    ===================================='''
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

'''
c1.Print("images/hmo/HMO_energy_dt_r{:.0f}_v{:.0f}_n10000000.pdf" % (tree4.vx0,tree4.vy0))
'''
