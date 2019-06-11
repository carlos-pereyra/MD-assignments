""" ============================================================
*  author: carlos p
*  purpose: particle in a box with open boundaries.
*
============================================================ """

from array import array
import numpy as np
from numpy import random
from datetime import datetime
random.seed(int(datetime.now().strftime("%s")))
import time
start_time=time.time()
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser

alpha_list  = [0.1,1,10]
temp_list   = [0.1,1,10]
L           = 20
natom       = 5
nstep       = int(10E3)
ndt         = 30
nal         = int(len(alpha_list))
ntemp       = int(len(temp_list))
#====================================
# tree1: TURFV
#====================================
t               	= array('f',[0])
ux              	= array('f',[0]*natom)
uy              	= array('f',[0]*natom)
x               	= array('f',[0]*natom)
y               	= array('f',[0]*natom)
fx              	= array('f',[0]*natom)
fy              	= array('f',[0]*natom)
vx              	= array('f',[2]*natom) #set initial velocity
vy              	= array('f',[3]*natom) #set initial velocity
# system
n               	= array('i',[0])
id              	= array('i',[0]*natom)
dt              	= array('f',[0])
dt_id           	= array('i',[0])
m               	= array('f',[1])
# energy per particle
ke_v            	= array('f',[0]*natom)
ke_u            	= array('f',[0]*natom)
u_lj            	= array('f',[0]*natom)
u_lj_pb         	= array('f',[0]*natom)
# gaussian number
alpha           	= array('f',[0])
alpha_id        	= array('i',[0])
xi_1            	= array('f',[0])
xi_2            	= array('f',[0])
eta_1           	= array('f',[0])
eta_2           	= array('f',[0])
beta_1         		= array('f',[0])
beta_2              = array('f',[0])
temp                = array('f',[0])
temp_id             = array('i',[0])
Kb                  = array('f',[1])
a			        = array('f',[0])
b			        = array('f',[0])
eku_true		    = array('f',[0])
#====================================
# tree2: Energy as a Function of time
#====================================
ke_v_natom_sum      = array('f',[0])
ke_u_natom_sum      = array('f',[0])
u_lj_natom_sum      = array('f',[0])
u_lj_pb_natom_sum   = array('f',[0])
ke_v_natom_avg      = array('f',[0])
ke_u_natom_avg      = array('f',[0])
u_lj_natom_avg      = array('f',[0])
u_lj_pb_natom_avg   = array('f',[0])
e_v_natom_avg       = array('f',[0]) #total average energy, velocity=v @t=n
e_u_natom_avg       = array('f',[0]) #total average energy, velocity=u @t=n
#====================================
# tree3: average energy over ntime
#====================================
ke_v_nsum		    = array('f',[0])
ke_u_nsum   	    = array('f',[0])
e_v_nsum    		= array('f',[0]) #average energy, velocity=v
e_u_nsum    		= array('f',[0]) #average energy, velocity=v
ke_v_navg   		= array('f',[0])
ke_u_navg   		= array('f',[0])
e_v_navg    		= array('f',[0]) #average energy, velocity=v
e_u_navg    		= array('f',[0]) #average energy, velocity=u
dt          		= array('f',[0])
#====================================
# tree4: System Parameters
#====================================
vx0                 = array('i',[int(vx[0]) ])
vy0             	= array('i',[int(vy[0]) ])
atoms           	= array('i',[natom])
n_time_steps    	= array('i',[0])
n_dt_steps      	= array('i',[0])

str0="data/dat_{}natom_{}nsteps_{:.0f}vx_{:.0f}vy_ndt_{}.root"
tfile=TFile(str0.format(natom,nstep,vx[0],vy[0],ndt),"recreate")
tree=TTree("timeSequence", "data with respect to time")
tree.Branch('atoms',atoms,'atoms/I')
tree.Branch('t',t,'t/F')
tree.Branch('ux',ux,'ux[atoms]/F')
tree.Branch('uy',uy,'uy[atoms]/F')
tree.Branch('x',x,'x[atoms]/F')
tree.Branch('y',y,'y[atoms]/F')
tree.Branch('fx',fx,'fx[atoms]/F')
tree.Branch('fy',fy,'fy[atoms]/F')
tree.Branch('vx',vx,'vx[atoms]/F')
tree.Branch('vy',vy,'vy[atoms]/F')
tree.Branch('id',id,'id[atoms]/I')
tree.Branch('ke_v',ke_v,'ke_v[atoms]/F')
tree.Branch('ke_u',ke_u,'ke_u[atoms]/F')
tree.Branch('u_lj',u_lj,'u_lj[atoms]/F')
tree.Branch('u_lj_pb',u_lj_pb,'u_lj_pb[atoms]/F')
tree.Branch('nstep',n,'n/I')
tree.Branch('dt',dt,'dt/F')
tree.Branch('dt_id',dt_id,'dt_id/I')
tree.Branch('alpha',alpha,'alpha/F')
tree.Branch('alpha_id',alpha_id,'alpha_id/I')
tree.Branch('temp',temp,'temp/F')
tree.Branch('temp_id'       ,temp_id,   'temp_id/I')
tree.Branch('m'             ,m,         'm/F')
#Tree 3
tree3=TTree("avgSystemEnergy", "energy vs dt")
tree3.Branch('ke_v_navg'    ,ke_v_navg, 'ke_v_navg/F')
tree3.Branch('ke_u_navg'    ,ke_u_navg, 'ke_u_navg/F')
tree3.Branch('e_v_navg'     ,e_v_navg,  'e_v_navg/F')
tree3.Branch('e_u_navg'     ,e_u_navg,  'e_u_navg/F')
tree3.Branch('eku_true'     ,eku_true,  'eku_true/F')
tree3.Branch('b'            ,b,         'b/F')
tree3.Branch('dt'           ,dt,        'dt/F')
tree3.Branch('dt_id'        ,dt_id,     'dt_id/I')
tree3.Branch('alpha'        ,alpha,     'alpha/F')
tree3.Branch('alpha_id'     ,alpha_id,  'alpha_id/I')
tree3.Branch('temp'         ,temp,      'temp/F')
tree3.Branch('temp_id'      ,temp_id,   'temp_id/I')
#Tree 4
tree4=TTree("systemParameters", "system parameters")
tree4.Branch('vx0'          ,vx0,       'vx0/I')
tree4.Branch('vy0'          ,vy0,       'vy0/I')

file = open('data/dat_{}natom_{}nsteps_{:.0f}vx_{:.0f}vy_ndt_{}'.format(natom,nstep,vx[0],vy[0],ndt),'w')
hdr="ek_v \t ek_u \t ep_v \t ep_u \t dt \t dt_id \t temp \t temp_id \t alpha \t alpha_id\n"
file.write(hdr)

def random_positions():
    #randomly generate x,y positions
    for i in range(0,len(x)):
        for j in range(i,len(x)):
            xij         = x[i]-x[j]
            yij         = y[i]-y[j]
            hyp_lj      = np.sqrt(xij**2 + yij**2)
            hyp_wall1   = np.sqrt((x[i]-L)**2+(y[i]-L)**2)
            hyp_wall2   = np.sqrt(x[i]**2+y[i]**2)
            while hyp_lj<1 or hyp_wall1<1 or hyp_wall2<1 or x[i]<1 or y[i]<1:
                x[i]=random.random()*(L-1)
                y[i]=random.random()*(L-1)
                xij         = x[i]-x[j]
                yij         = y[i]-y[j]
                hyp_lj      = np.sqrt(xij**2 + yij**2)
                hyp_wall1   = np.sqrt((x[i]-L)**2+(y[i]-L)**2)
                hyp_wall2   = np.sqrt(x[i]**2+y[i]**2)
                if(i==j):
                    break


def wall_potential(xi,nu):
    #wall scalar potential
    u_left_boundary     = ((1/xi)**12-2*(1/xi)**6)
    u_right_boundary    = ((1/(L-xi))**12-2*(1/(L-xi))**6)
    u_bottom_boundary   = ((1/nu)**12-2*(1/nu)**6)
    u_top_boundary      = ((1/(L-nu))**12-2*(1/(L-nu))**6)
    u=u_left_boundary+u_right_boundary+u_bottom_boundary+u_top_boundary
    return u

def wall_force(xi,nu):
    #wall force vector force
    fx_left_boundary    = 12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary   = -12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx                  = fx_left_boundary+fx_right_boundary
    
    fy_bottom_boundary  = 12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary     = -12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy                  = fy_bottom_boundary+fy_top_boundary
    return fx,fy

def lj_potential(xij,yij,i=0,j=0,n=0):
    #lennard jones scalar potential
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2-body separation distance
    if r<ri:
        if r==0:
            label="lj_potential()->avoid divide by zero [i={} j={} n={}]"
            print(label.format(i,j,n))
            u = 10000
        else:
            u = (1/r)**12-2*(1/r)**6
    if ri<=r<rc:
        u = -6.129*(r-rc)**2 + -4.655*(r-rc)**3
    elif r>=rc:
        u = 0

    return u

'''
def lj_force(xij,yij):
    #lennard jones force vector
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    if r<ri:
        if r==0:
            print("lj_force()->avoid divide by zero")
            fx = 10000
            fy = 10000
        else:
            fx = (1/r**13 + 1/r**7)*(xij/r) #*(12/np.sqrt(2))
            fy = (1/r**13 + 1/r**7)*(xij/r) #*(12/np.sqrt(2))
    if ri<=r<rc:
        fx      = 12.258*(xij-rc)**1 + 13.966*(xij-rc)**2
        fy      = 12.258*(yij-rc)**1 + 13.966*(yij-rc)**2
    elif r>=rc:
        fx      = 0
        fy      = 0
    
    return fx,fy'''


def lj_force(xij,yij):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    if r is 0:
        r+=.0001
    ir = 1/r
    
    if r<ri:
        fx = (12/r**13 - 12/r**7)*(xij*ir)
        fy = (12/r**13 - 12/r**7)*(xij*ir)
        lj = 1
        sp = 0
    
    if ri<=r<rc:
        fx = 12.258*((r-rc)**1)*(xij*ir) + 13.966*((r-rc)**2)*(xij*ir)
        fy = 12.258*((r-rc)**1)*(yij*ir) + 13.966*((r-rc)**2)*(yij*ir)
        lj = 0
        sp = 1
    elif r>=rc:
        fx = 0
        fy = 0
        lj = 0
        sp = 0
    
    return fx, fy


def main():
    random_positions()
    random_positions()
    print("randomization done --- %s seconds --- changing" % (time.time()-start_time))

    #initialize matrices
    lj_u_mat=np.zeros((natom,natom))
    lj_pb_u_mat=np.zeros((natom,natom))
    lj_f_mat=np.zeros((2,natom,natom)) #mat[i][j][0] fx, mat[i][j][0] fy
    lj_pb_f_mat=np.zeros((2,natom,natom)) #periodic boundary forces
    for n_al in range(0,nal):
        #============================
        #   Alpha coef Loop
        #============================
        alpha[0]    = alpha_list[n_al]
        alpha_id[0] = n_al
        for n_temp in range(0,ntemp):
            #============================
            #   Temperature Loop
            #============================
            temp[0]     = temp_list[n_temp]
            temp_id[0]  = n_temp
        
            for n_dt in range(0,ndt):
                #============================
                #   dt Loop
                #============================
                random_positions()
                random_positions()
                dt[0]       = 0.001*(n_dt+1)
                dt_id[0]    = n_dt
                ke_v_nsum   = array('f',[0])
                ke_u_nsum   = array('f',[0])
                e_v_nsum    = array('f',[0])
                e_u_nsum    = array('f',[0])
                
                inv2m       = 1/(2.*m[0])
                a[0]        = (1-alpha[0]*dt[0]*inv2m )/(1+alpha[0]*dt[0]*inv2m )
                b[0]        = 1/(1+alpha[0]*dt[0]*inv2m )
                c           = a[0]/np.sqrt(b[0])
                d           = np.sqrt(b[0])
                for n_time in range(0,nstep):
                    #============================
                    #  Time Loop
                    #============================
                    ke_v_natom_sum      = array('f',[0])
                    ke_u_natom_sum      = array('f',[0])
                    u_lj_natom_sum      = array('f',[0])
                    u_lj_pb_natom_sum   = array('f',[0])
                    
                    #calculate force on each particle
                    for i in range(0,natom):
                        #ith particle
                        for j in range(i,natom):
                            #============================
                            #  Lennard-Jones Loop
                            #============================
                            xij=x[i]-x[j]
                            yij=y[i]-y[j]
                            
                            #interbox collision
                            if i==j: #diagonal
                                lj_u_mat[i][j]      = 0 # u
                                lj_f_mat[0][i][j]   = 0 # fx
                                lj_f_mat[1][i][j]   = 0 # fy
                            else: #off diagonal
                                lj_u_mat[i][j]      = lj_potential(xij,yij,i,j,n_time)
                                lj_u_mat[j][i]      = lj_u_mat[i][j]
                                f_x,f_y             = lj_force(xij,yij)
                                lj_f_mat[0][i][j]   = f_x
                                lj_f_mat[1][i][j]   = f_y
                                lj_f_mat[0][j][i]   =-f_x
                                lj_f_mat[1][j][i]   =-f_y
                            
                            # periodic boundary conditions
                            if xij>0.5*L:   xij-=L
                            if xij<=-0.5*L: xij+=L
                            if yij>0.5*L:   yij-=L
                            if yij<=-0.5*L: yij+=L
                            if i==j: #diagonal
                                lj_pb_u_mat[i][j]    = 0 # u
                                lj_pb_f_mat[0][i][j] = 0 # fx
                                lj_pb_f_mat[1][i][j] = 0 # fy
                            else: #off diagonal
                                lj_pb_u_mat[i][j]    = lj_potential(xij,yij)
                                lj_pb_u_mat[j][i]    = lj_pb_u_mat[i][j]
                                f_x,f_y              = lj_force(xij,yij)
                                lj_pb_f_mat[0][i][j] = f_x
                                lj_pb_f_mat[1][i][j] = f_y
                                lj_pb_f_mat[0][j][i] =-f_x
                                lj_pb_f_mat[1][j][i] =-f_y
                        
                    #update Positions
                    for i in range(0,natom):
                        # Langevin Motion Algorithm
                        xi_1[0]=random.uniform(0,1) # xi1 calculation
                        if xi_1[0]==0:
                            print("change xi-1: {}".format(xi_1[0]))
                            xi_1[0]=1-xi_1[0]
                        
                        xi_2[0]=random.uniform(0,1) # xi2 calculation
                        if xi_2[0]==0:
                            print("change xi-2: {}".format(xi_2[0]))
                            xi_2[0]=1-xi2[0]
                        
                        # gaussian number process
                        eta_1[0]=np.sqrt(-2*np.log(xi_1[0]))*np.cos(2*np.pi*xi_2[0])
                        eta_2[0]=np.sqrt(-2*np.log(xi_1[0]))*np.sin(2*np.pi*xi_2[0])
                        beta_1[0]=np.sqrt(2*alpha[0]*Kb[0]*temp[0]*dt[0])*eta_1[0]
                        beta_2[0]=np.sqrt(2*alpha[0]*Kb[0]*temp[0]*dt[0])*eta_2[0]
                        # turfv
                        t[0]=dt[0]*(n_time+1)
                        ux[i]   = d*(vx[0]+(dt[0]*fx[i]+beta_1[0])*inv2m )
                        uy[i]   = d*(vy[0]+(dt[0]*fy[i]+beta_2[0])*inv2m )
                        x[i]    = x[i]+d*dt[0]*ux[i]
                        y[i]    = y[i]+d*dt[0]*uy[i]
                            
                        #reenter box boundary
                        if x[i]>=L:     x[i]-=L
                        if x[i]<0:      x[i]+=L
                        if y[i]>=L:     y[i]-=L
                        if y[i]<0:      y[i]+=L
                        
                        f_lj_x=np.sum(lj_f_mat[0][i])
                        f_lj_y=np.sum(lj_f_mat[1][i])
                        f_pb_x=np.sum(lj_pb_f_mat[0][i])
                        f_pb_y=np.sum(lj_pb_f_mat[1][i])
                        fx[i]=f_lj_x+f_pb_x
                        fy[i]=f_lj_y+f_pb_y
                        #vx[i]=ux[i]+dt[0]*fx[i]/(2*m[0])
                        #vy[i]=uy[i]+dt[0]*fy[i]/(2*m[0])
                        vx[i]=c*ux[0]+(dt[0]*fx[i]+beta_1[0])*inv2m
                        vy[i]=c*uy[0]+(dt[0]*fy[i]+beta_2[0])*inv2m

                        # i-th Particle Energy Calculation
                        ke_v[i]     = 0.5*m[0]*(vx[i]**2+vy[i]**2)
                        ke_u[i]     = 0.5*m[0]*(ux[i]**2+uy[i]**2)
                        u_lj[i]     = np.sum(lj_u_mat[i])
                        u_lj_pb[i]  = np.sum(lj_pb_u_mat[i])
                        id[i]=i
                        
                        # sum paritlce energy calculation
                        ke_v_natom_sum[0]   += ke_v[i] #eKn_sum_v
                        ke_u_natom_sum[0]   += ke_u[i]
                        u_lj_natom_sum[0]   += u_lj[i]
                        u_lj_pb_natom_sum[0]+= u_lj_pb[i]
                        itrstr="i-> {} ux->{} uy->{} vx->{} vy->{} fx->{} fy->{} beta1->{} eta->{} alpha->{} temp->{} dt->{}"
                        #print(itrstr.format(i,ux[i],uy[i],vx[i],vy[i],fx[i],fy[i],beta_1[0],eta_1[0],alpha[0],temp[0],dt[0]))
                        #if ux[i]>10000: exit(1)
                    
                    # average energy per particle
                    u_lj_natom_avg[0]       = u_lj_natom_sum[0]/float(natom)
                    u_lj_pb_natom_avg[0]    = u_lj_pb_natom_sum[0]/float(natom)
                    ke_v_natom_avg[0]       = ke_v_natom_sum[0]/float(natom)
                    ke_u_natom_avg[0]       = ke_u_natom_sum[0]/float(natom)
                    e_v_natom_avg[0]        = ke_v_natom_avg[0]+u_lj_natom_avg[0]+u_lj_pb_natom_avg[0]
                    e_u_natom_avg[0]        = ke_u_natom_avg[0]+u_lj_natom_avg[0]+u_lj_pb_natom_avg[0]

                    # Sum averaged energies
                    ke_v_nsum[0]  += ke_v_natom_avg[0]
                    ke_u_nsum[0]  += ke_u_natom_avg[0]
                    e_v_nsum[0]   += e_v_natom_avg[0]
                    e_u_nsum[0]   += e_u_natom_avg[0]
            
                    if n_time<1000:
                        tree.Fill()

                # average energy per time
                ke_v_navg[0] = ke_v_nsum[0]/float(nstep)    #<KE_p^n>
                ke_u_navg[0] = ke_u_nsum[0]/float(nstep)    #<KE_p^n>
                e_v_navg[0]  = e_v_nsum[0]/float(nstep)     #<E_p^n>
                e_u_navg[0]  = e_u_nsum[0]/float(nstep)     #<E_p^n>
                eku_true[0]  = 0.5*Kb[0]*temp[0]*b[0]
                timestr      = time.time()-start_time
                tree3.Fill() #energy vs dt (at alpha=1 and temp=1)
                str1="{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.3f}\t{:.0f}\t{:.1f}\t{:.0f}\t{:.1f}\t{:.0f}\n"
                file.write(str1.format(ke_v_navg[0],
                                       ke_u_navg[0],
                                       e_v_navg[0],
                                       e_u_navg[0],
                                       dt[0],
                                       dt_id[0],
                                       temp[0],
                                       temp_id[0],
                                       alpha[0],
                                       alpha_id[0]))
                          
                print("=============================================")
                print("ndt: %s -- %s seconds -- %s alpha -- %s temp" % (dt[0],timestr,alpha[0],temp[0]))
                str2="ekv->{:.2f} \t eku->{:.2f} \t epv->{:.2f} \t epu->{:.2f} \t dt->{:.3f}\t dtId->{:.0f} \t temp->{:.1f} \t tempId->{:.0f} \t alpha->{:.1f}\t alpha_id->{:.0f}\n"
                print(str2.format(ke_v_navg[0],
                                  ke_u_navg[0],
                                  e_v_navg[0],
                                  e_u_navg[0],
                                  dt[0],
                                  dt_id[0],
                                  temp[0],
                                  temp_id[0],
                                  alpha[0],
                                  alpha_id[0]))
                    
    tree4.Fill()
    # Write2Disk
    tree.Write() #data respect to time
    tree3.Write() # energy vs dt
    tree4.Write() # system parameters
    tfile.Close()
    file.close()

if __name__ == "__main__":
    main()
