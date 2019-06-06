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
import fileOut as f
L=f.l

def random_positions():
    for i in range(0,len(f.x)):
        for j in range(i,len(f.x)):
            while (np.sqrt((f.x[i]-f.x[j])**2+(f.y[i]-f.y[j])**2)<1) or (np.sqrt((f.x[i]-f.l)**2+(f.y[i]-f.l)**2)<1) or (np.sqrt(f.x[i]**2+f.y[i]**2)<1) or f.x[i]<1 or f.y[i]<1:
                #print("too close together i {} j {}".format(i,j))
                f.x[i]=random.random()*(f.l-1)
                f.y[i]=random.random()*(f.l-1)
                if(i==j):
                    break


def wall_potential(xi,nu):
    u_left_boundary     = ((1/xi)**12-2*(1/xi)**6)
    u_right_boundary    = ((1/(L-xi))**12-2*(1/(L-xi))**6)
    u_bottom_boundary   = ((1/nu)**12-2*(1/nu)**6)
    u_top_boundary      = ((1/(L-nu))**12-2*(1/(L-nu))**6)
    u=u_left_boundary+u_right_boundary+u_bottom_boundary+u_top_boundary
    return u


def lennard_jones_potential(r):
    u_lj=(1/r)**12-2*(1/r)**6
    return u_lj


def wall_force(xi,nu):
    fx_left_boundary    = 12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary   = -12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx                  = fx_left_boundary+fx_right_boundary
    
    fy_bottom_boundary  = 12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary     = -12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy                  = fy_bottom_boundary+fy_top_boundary
    return fx,fy


def lennard_jones_force(xij,yij):
    gamma = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
    fx    = (1/gamma**13 + 1/gamma**7)*(xij/gamma)*(12/np.sqrt(2))
    fy    = (1/gamma**13 + 1/gamma**7)*(yij/gamma)*(12/np.sqrt(2))
    return fx,fy


def main():
    # random positioning
    random_positions()
    random_positions()
    print("randomization done --- %s seconds --- changing" % (time.time()-start_time))

    #initialize matrices
    r=np.zeros((f.natom,f.natom))
    u_jones_matrix=np.zeros((f.natom,f.natom))
    f_jones_matrix_x=np.zeros((f.natom,f.natom))
    f_jones_matrix_y=np.zeros((f.natom,f.natom))

    for n_al in range(0,f.nal): #============================
                                #   Alpha coef Loop
                                #============================
        f.alpha[0]    = f.alpha_list[n_al]
        f.alpha_id[0] = n_al
        
        for n_temp in range(0,f.ntemp): #============================
                                        #   Temperature Loop
                                        #============================
            f.temp[0]     = f.temp_list[n_temp]
            f.temp_id[0]  = n_temp
        
            for n_dt in range(0,f.ndt): #============================
                                        #   dt Loop
                                        #============================
                random_positions()
                random_positions()
                f.dt[0]             = 0.001*(n_dt+1)
                f.dt_id[0]          = n_dt
                f.total_ke_v_sum    = array('f',[0])
                f.total_ke_u_sum    = array('f',[0])
                f.total_v_sum       = array('f',[0])
                f.total_u_sum       = array('f',[0])
                print("ndt: %s -- %s seconds -- %s alpha -- %s temp" % (f.dt[0],
                                                               time.time()-start_time,
                                                               f.alpha[0],
                                                               f.temp[0]))
                
                for n_time in range(0,f.nstep): #============================
                                                #  Time Loop
                                                #============================
                    f.ke_v_sum_i    = array('f',[0])
                    f.ke_u_sum_i    = array('f',[0])
                    f.u_jones_sum_i = array('f',[0])
                    f.u_wall_sum_i  = array('f',[0])
                    
                    for i in range(0,f.natom): #ith particle
                        
                        for j in range(i,f.natom): #============================
                                                   #  Lennard-Jones Loop
                                                   #============================
                            xij=f.x[i]-f.x[j]
                            yij=f.y[i]-f.y[j]
                            r[i][j]=(xij**2+yij**2)**0.5
                            r[j][i]=r[i][j]
                            if i==j:
                                r[i][j]=0
                                u_jones_matrix[i][j]=0
                                f_jones_matrix_x[i][j]=0
                                f_jones_matrix_y[i][j]=0
                            else:
                                u_jones_matrix[i][j]=lennard_jones_potential(r[i][j])
                                u_jones_matrix[j][i]=u_jones_matrix[i][j]
                                f_jones_matrix_x[i][j],f_jones_matrix_y[i][j]=lennard_jones_force(xij,yij)
                                f_jones_matrix_x[j][i]=-f_jones_matrix_x[i][j]
                                f_jones_matrix_y[j][i]=-f_jones_matrix_y[i][j]
                    
                        '''========================================
                            Verlet Motion Algorithm
                            ========================================'''
                        '''
                        f.t[0]=f.dt[0]*n_time
                        f.ux[i]=f.vx[i]+f.dt[0]*f.fx[i]/(2*f.m[0])
                        f.uy[i]=f.vy[i]+f.dt[0]*f.fy[i]/(2*f.m[0])
                        f.x[i]=f.x[i]+f.ux[i]*f.dt[0]
                        f.y[i]=f.y[i]+f.uy[i]*f.dt[0]'''
                        
                        '''========================================
                            Langevin Motion Algorithm
                            ========================================'''
                        f.xi_1[0]=random.uniform(0,1) # xi1 calculation
                        if f.xi_1[0]==0:
                            print("change xi-1: {}".format(f.xi_1[0]))
                            f.xi_1[0]=1-f.xi_1[0]
                        
                        f.xi_2[0]=random.uniform(0,1) # xi2 calculation
                        if f.xi_2[0]==0:
                            print("change xi-2: {}".format(f.xi_2[0]))
                            f.xi_2[0]=1-f.xi2[0]
                        
                        # random number process
                        f.eta_1[0]=np.sqrt(-2*np.log(f.xi_1[0]))*np.cos(2*np.pi*f.xi_2[0])
                        f.eta_2[0]=np.sqrt(-2*np.log(f.xi_1[0]))*np.sin(2*np.pi*f.xi_2[0])
                        f.beta_1[0]=np.sqrt(2*f.alpha[0]*f.Kb[0]*f.temp[0]*f.dt[0])*f.eta_1[0]
                        f.beta_2[0]=np.sqrt(2*f.alpha[0]*f.Kb[0]*f.temp[0]*f.dt[0])*f.eta_2[0]
                        a=(1-f.alpha[0]*f.dt[0]/(2*f.m[0]) )/(1+f.alpha[0]*f.dt[0]/(2*f.m[0]) )
                        b=1/(1+f.alpha[0]*f.dt[0]/(2*f.m[0]) )
                        # turfv
                        f.t[0]=f.dt[0]*(n_time+1)
                        f.ux[i]=np.sqrt(b)*(f.vx[0]+(f.dt[0]*f.fx[i]+f.beta_1[0])/(2*f.m[0]) )
                        f.uy[i]=np.sqrt(b)*(f.vy[0]+(f.dt[0]*f.fy[i]+f.beta_2[0])/(2*f.m[0]) )
                        f.x[i]=f.x[i]+np.sqrt(b)*f.dt[0]*f.ux[i]
                        f.y[i]=f.y[i]+np.sqrt(b)*f.dt[0]*f.uy[i]

                        '''================
                            lennard jones
                            forces
                            ================'''
                        f_lj_x=np.sum(f_jones_matrix_x[i])
                        f_lj_y=np.sum(f_jones_matrix_y[i])
                        #f_wall_x,f_wall_y=wall_force(f.x[i],f.y[i])
                        '''================
                            periodic boundary
                            forces
                            ================'''
                        xi=f.x[i]
                        xj=-f.x[i]
                        xij=xj-xi
                        if xij>=0.5*L: xij-=L
                        if xij<-0.5*L: xij+=L
                        
                        yi=f.y[i]
                        yj=-f.y[i]
                        yij=yj-yi
                        if yij>=0.5*L: yij-=L
                        if yij<-0.5*L: yij+=L
                        
                        f_pdc_x,f_pdc_y=lennard_jones_force(xij,yij)
                        
                        '''================
                            total forces
                            ================'''
                        f.fx[i]=f_pdc_x+f_lj_x
                        f.fy[i]=f_pdc_y+f_lj_y
                        #f.fx[i]=f_wall_x+f_lj_x
                        #f.fy[i]=f_wall_y+f_lj_y
                        
                        '''================
                            update velocity
                            ================'''
                        f.vx[i]=f.ux[i]+f.dt[0]*f.fx[i]/(2*f.m[0])
                        f.vy[i]=f.uy[i]+f.dt[0]*f.fy[i]/(2*f.m[0])
                        
                        #========================================
                        # i-th Particle Energy Calculation
                        #========================================
                        f.ke_v[i]=0.5*f.m[0]*(f.vx[i]**2+f.vy[i]**2)
                        f.ke_u[i]=0.5*f.m[0]*(f.ux[i]**2+f.uy[i]**2)
                        f.u_jones[i]=np.sum(u_jones_matrix[i])
                        f.u_wall[i]=wall_potential(f.x[i],f.y[i])
                        f.id[i]=i
                        #========================================
                        # SUM Paritlce Energy Calculation
                        #========================================
                        f.ke_v_sum_i[0]     += f.ke_v[i]
                        f.ke_u_sum_i[0]     += f.ke_u[i]
                        f.u_jones_sum_i[0]  += f.u_jones[i]
                        f.u_wall_sum_i[0]   += f.u_wall[i]
                    
                    #========================================
                    # TimeLoop: Average Paritlce Energy Calculation
                    #========================================
                    f.ke_v_avg_i[0]     = f.ke_v_sum_i[0]/float(f.natom)
                    f.ke_u_avg_i[0]     = f.ke_u_sum_i[0]/float(f.natom)
                    f.u_jones_avg_i[0]  = f.u_jones_sum_i[0]/(f.natom)
                    f.u_wall_avg_i[0]   = f.u_wall_sum_i[0]/(f.natom)
                    f.total_v_i[0]      = f.ke_v_avg_i[0]+f.u_wall_avg_i[0]+f.u_jones_avg_i[0]
                    f.total_u_i[0]      = f.ke_u_avg_i[0]+f.u_wall_avg_i[0]+f.u_jones_avg_i[0]
                    
                    #print("ke u: {:.2f} ke v: {:.2f} keU sum: {:.2f} keV sum: {:.2f} {} {}".format(f.ke_u_avg_i[0],f.ke_u_avg_i[0],f.ke_u_sum_i[0],f.ke_v_sum_i[0],f.natom, f.nstep))
                    #========================================
                    # Triggering condition
                    #========================================
                    if n_time<f.nrecord:
                        f.tree.Fill()
                        f.tree2.Fill() #energy per time

                    #========================================
                    # Prep for energy vs dt
                    # - sum energies
                    #========================================
                    f.total_ke_v_sum[0] += f.ke_v_avg_i[0]
                    f.total_ke_u_sum[0] += f.ke_u_avg_i[0]
                    f.total_v_sum[0]    += f.total_v_i[0] #avg particle energy sum, velocity=v
                    f.total_u_sum[0]    += f.total_u_i[0] #avg particle energy sum, velocity=u

                #========================================
                # dt Loop averageEnergy(dt)
                #========================================
                f.total_ke_v_avg[0] = f.total_ke_v_sum[0]/float(f.nstep)    #<KE_p^n>
                f.total_ke_u_avg[0] = f.total_ke_u_sum[0]/float(f.nstep)    #<KE_p^n>
                f.total_e_v_avg[0]  = f.total_v_sum[0]/float(f.nstep)       #<E_p^n>
                f.total_e_u_avg[0]  = f.total_u_sum[0]/float(f.nstep)       #<E_p^n>
                f.tree3.Fill() #energy vs dt (at alpha=1 and temp=1)

                #print("AllTime SUM -> keU: {:.2f} keV: {:.2f} totalU: {:.2f} keV sum: {:.2f} {} {}".format(f.total_ke_u_sum[0],f.total_ke_u_sum[0],f.total_ke_u_sum[0],f.total_ke_v_sum[0],f.natom, f.nstep))
                
                print("AllTime AVG -> keU: {:.2f} keV: {:.2f} totalU: {:.2f} totalV: {:.2f} {} {}".format(f.total_ke_u_avg[0],f.total_ke_v_avg[0],f.total_e_u_avg[0],f.total_e_v_avg[0],f.natom, f.nstep))

            f.atoms[0]              = f.natom
            f.n_time_steps[0]       = f.nstep
            f.n_dt_steps[0]         = f.ndt
            f.tree4.Fill()

    '''======================== Write2Disk ==========================='''
    f.tree.Write() # everything respect to time
    f.tree2.Write() # energy vs time
    f.tree3.Write() # energy vs dt
    f.tree4.Write() # system parameters
    f.file.Close()

if __name__ == "__main__":
    main()
