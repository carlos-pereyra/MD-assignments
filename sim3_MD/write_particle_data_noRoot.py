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

alpha_list  = [0,0.1,1,10]
temp_list   = [0.1,1,10]
L           = 20
natom       = 5
nstep       = int(10E0)
ndt         = 20
nal         = int(len(alpha_list))
ntemp       = int(len(temp_list))
#====================================
# tree1: TURFV
#====================================
t               = array('f',[0])
ux              = array('f',[0]*natom)
uy              = array('f',[0]*natom)
x               = array('f',[0]*natom)
y               = array('f',[0]*natom)
fx              = array('f',[0]*natom)
fy              = array('f',[0]*natom)
vx              = array('f',[2]*natom) #set initial velocity
vy              = array('f',[3]*natom) #set initial velocity
# system
n               = array('i',[0])
id              = array('i',[0]*natom)
dt              = array('f',[0])
dt_id           = array('i',[0])
m               = array('f',[1])
# energy per particle
ke_v            = array('f',[0]*natom)
ke_u            = array('f',[0]*natom)
u_wall          = array('f',[0]*natom)
u_jones         = array('f',[0]*natom)
# gaussian number
alpha           = array('f',[0])
alpha_id        = array('i',[0])
xi_1            = array('f',[0])
xi_2            = array('f',[0])
eta_1           = array('f',[0])
eta_2           = array('f',[0])
beta_1          = array('f',[0])
beta_2          = array('f',[0])
temp            = array('f',[0])
temp_id         = array('i',[0])
Kb              = array('f',[1])
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
ke_v_nsum   = array('f',[0])
ke_u_nsum   = array('f',[0])
e_v_nsum    = array('f',[0]) #average energy, velocity=v
e_u_nsum    = array('f',[0]) #average energy, velocity=v
ke_v_navg   = array('f',[0])
ke_u_navg   = array('f',[0])
e_v_navg    = array('f',[0]) #average energy, velocity=v
e_u_navg    = array('f',[0]) #average energy, velocity=u
dt          = array('f',[0])
#====================================
# tree4: System Parameters
#====================================
vx0             =   array('i',[int(vx[0]) ])
vy0             =   array('i',[int(vy[0]) ])
atoms           =   array('i',[natom])
n_time_steps    =   array('i',[0])
n_dt_steps      =   array('i',[0])

file = open('data/dat_{}natom_{}nsteps_{:.0f}vx_{:.0f}vy'.format(natom,nstep,vx[0],vy[0]),'w')
file.write("ek_v\tek_u\tep_v\tep_u\tdt\ttemp\talpha\tk\n")

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

def lj_potential(xij,yij):
    #lennard jones scalar potential
    hyp = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    u_lj=(1/hyp)**12-2*(1/hyp)**6
    return u_lj

def lj_force(xij,yij):
    #lennard jones force vector
    hyp = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    fx  = (1/hyp**13 + 1/hyp**7)*(xij/hyp)*(12/np.sqrt(2))
    fy  = (1/hyp**13 + 1/hyp**7)*(yij/hyp)*(12/np.sqrt(2))
    return fx,fy


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
                
                for n_time in range(0,nstep):
                    #============================
                    #  Time Loop
                    #============================
                    ke_v_natom_sum      = array('f',[0])
                    ke_u_natom_sum      = array('f',[0])
                    u_lj_natom_sum      = array('f',[0])
                    u_lj_pb_natom_sum   = array('f',[0])
                    
                    for i in range(0,natom):
                        #ith particle
                        
                        for j in range(i,natom):
                            #============================
                            #  Lennard-Jones Loop
                            #============================
                            xij=x[i]-x[j]
                            yij=y[i]-y[j]
                            #inter-box collision
                            if i==j:
                                lj_u_mat[i][j]      = 0 # u
                                lj_f_mat[0][i][j]   = 0 # fx
                                lj_f_mat[1][i][j]   = 0 # fy
                            else:
                                lj_u_mat[i][j]      = lj_potential(xij,yij)
                                lj_u_mat[j][i]      = lj_u_mat[i][j]
                                f_x,f_y             = lj_force(xij,yij)
                                lj_f_mat[0][i][j]   = f_x
                                lj_f_mat[1][i][j]   = f_y
                                lj_f_mat[0][j][i]   =-f_x
                                lj_f_mat[1][j][i]   =-f_y
                            
                            # periodic boundary jth-atom mirror image
                            if xij>0.5*L:   xij-=L
                            if xij<=-0.5*L: xij+=L
                            if yij>0.5*L:   yij-=L
                            if yij<=-0.5*L: yij+=L
                            #periodic boundary forces & potential
                            if i==j:
                                lj_pb_u_mat[i][j]    = 0 # u
                                lj_pb_f_mat[0][i][j] = 0 # fx
                                lj_pb_f_mat[1][i][j] = 0 # fy
                            else:
                                lj_pb_u_mat[i][j]    = lj_potential(xij,yij)
                                lj_pb_u_mat[j][i]    = lj_pb_u_mat[i][j]
                                f_x,f_y              = lj_force(xij,yij)
                                lj_pb_f_mat[0][i][j] = f_x
                                lj_pb_f_mat[1][i][j] = f_y
                                lj_pb_f_mat[0][j][i] =-f_x
                                lj_pb_f_mat[1][j][i] =-f_y
                            #reenter box boundary
                            if x[j]>=L:     x[j]-=L
                            if x[j]<0:      x[j]+=L
                            if y[j]>=L:     y[j]-=L
                            if y[j]<0:      y[j]+=L
                        
                        '''# ========================================
                           # Langevin Motion Algorithm
                           #
                           # ========================================'''
                        xi_1[0]=random.uniform(0,1) # xi1 calculation
                        if xi_1[0]==0:
                            print("change xi-1: {}".format(xi_1[0]))
                            xi_1[0]=1-xi_1[0]
                        
                        xi_2[0]=random.uniform(0,1) # xi2 calculation
                        if xi_2[0]==0:
                            print("change xi-2: {}".format(xi_2[0]))
                            xi_2[0]=1-xi2[0]
                        
                        # random number process
                        eta_1[0]=np.sqrt(-2*np.log(xi_1[0]))*np.cos(2*np.pi*xi_2[0])
                        eta_2[0]=np.sqrt(-2*np.log(xi_1[0]))*np.sin(2*np.pi*xi_2[0])
                        beta_1[0]=np.sqrt(2*alpha[0]*Kb[0]*temp[0]*dt[0])*eta_1[0]
                        beta_2[0]=np.sqrt(2*alpha[0]*Kb[0]*temp[0]*dt[0])*eta_2[0]
                        a=(1-alpha[0]*dt[0]/(2*m[0]) )/(1+alpha[0]*dt[0]/(2*m[0]) )
                        b=1/(1+alpha[0]*dt[0]/(2*m[0]) )
                        # turfv
                        t[0]=dt[0]*(n_time+1)
                        ux[i]=np.sqrt(b)*(vx[0]+(dt[0]*fx[i]+beta_1[0])/(2*m[0]) )
                        uy[i]=np.sqrt(b)*(vy[0]+(dt[0]*fy[i]+beta_2[0])/(2*m[0]) )
                        x[i]=x[i]+np.sqrt(b)*dt[0]*ux[i]
                        y[i]=y[i]+np.sqrt(b)*dt[0]*uy[i]
                        f_lj_x=np.sum(lj_f_mat[0][i])
                        f_lj_y=np.sum(lj_f_mat[1][i])
                        f_pb_x=np.sum(lj_pb_f_mat[0][i])
                        f_pb_y=np.sum(lj_pb_f_mat[1][i])
                        fx[i]=f_pb_x+f_lj_x
                        fy[i]=f_pb_y+f_lj_y
                        vx[i]=ux[i]+dt[0]*fx[i]/(2*m[0])
                        vy[i]=uy[i]+dt[0]*fy[i]/(2*m[0])
                        
                        # i-th Particle Energy Calculation
                        ke_v[i]         = 0.5*m[0]*(vx[i]**2+vy[i]**2)
                        ke_u[i]         = 0.5*m[0]*(ux[i]**2+uy[i]**2)
                        u_jones[i]      = np.sum(lj_u_mat[i])
                        u_periodic      = np.sum(lj_pb_u_mat[i])
                        id[i]=i
                        
                        # sum paritlce energy calculation
                        ke_v_natom_sum[0]   += ke_v[i] #eKn_sum_v
                        ke_u_natom_sum[0]   += ke_u[i]
                        u_lj_natom_sum[0]   += u_jones[i]
                        u_lj_pb_natom_sum[0]+= u_periodic
                    
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

                # average energy per time
                ke_v_navg[0] = ke_v_nsum[0]/float(nstep)    #<KE_p^n>
                ke_u_navg[0] = ke_u_nsum[0]/float(nstep)    #<KE_p^n>
                e_v_navg[0]  = e_v_nsum[0]/float(nstep)     #<E_p^n>
                e_u_navg[0]  = e_u_nsum[0]/float(nstep)     #<E_p^n>
                file.write("{:.2f}\t{:.2f}\t{:.2f}\t{:.2f}\t{:.3f}\t{:.1f}\t{:.1f}\n".format(ke_v_navg[0],
                                                                                           ke_u_navg[0],
                                                                                           e_v_navg[0],
                                                                                           e_u_navg[0],
                                                                                           dt[0],
                                                                                           temp[0],
                                                                                           alpha[0]))
                          
                timestr=time.time()-start_time
                print("=============================================")
                print("ndt: %s -- %s seconds -- %s alpha -- %s temp" % (dt[0],timestr,alpha[0],temp[0]))
                print("ek_u: {:.2f} ek_v: {:.2f} e_v: {:.2f} e_u: {:.2f} {} {} {}".format(ke_v_navg[0],
                                                                                          ke_u_navg[0],
                                                                                          e_v_navg[0],
                                                                                          e_u_navg[0],
                                                                                          dt[0],
                                                                                          temp[0],alpha[0]))

    # Write2Disk
    file.close()

if __name__ == "__main__":
    main()
