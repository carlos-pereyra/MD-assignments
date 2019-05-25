""" ============================================================
*  author: carlos p
*  purpose: particle in a box
*
============================================================ """

from array import array
import numpy as np
from numpy import random
from datetime import datetime
import time
start_time=time.time()
import fileOut as f
l=f.l

random.seed(int(datetime.now().strftime("%s")))

def random_positions():
    for i in range(0,len(f.x_dum)):
        for j in range(i,len(f.x_dum)):
            while (np.sqrt((f.x_dum[i]-f.x_dum[j])**2+(f.y_dum[i]-f.y_dum[j])**2)<1) or (np.sqrt((f.x_dum[i]-f.l)**2+(f.y_dum[i]-f.l)**2)<1) or (np.sqrt(f.x_dum[i]**2+f.y_dum[i]**2)<1) or f.x_dum[i]<1 or f.y_dum[i]<1:
                #print("too close together i {} j {}".format(i,j))
                f.x_dum[i]=random.random()*(f.l-1)
                f.y_dum[i]=random.random()*(f.l-1)
                if(i==j):
                    break

    f.x_dum[0]=1
    f.y_dum[0]=1


def wall_potential(xi,nu):
    u_left_boundary     = ((1/xi)**12-2*(1/xi)**6)
    u_right_boundary    = ((1/(l-xi))**12-2*(1/(l-xi))**6)
    u_bottom_boundary   = ((1/nu)**12-2*(1/nu)**6)
    u_top_boundary      = ((1/(l-nu))**12-2*(1/(l-nu))**6)
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


def interaction_force(xi,nu):
    fx_left_boundary    = 12*((1/xi)**13-(1/xi)**7)
    fx_right_boundary   = -12*((1/(l-xi))**13-(1/(l-xi))**7)
    fx                  = fx_left_boundary+fx_right_boundary
    
    fy_bottom_boundary  = 12*((1/nu)**13-(1/nu)**7)
    fy_top_boundary     = -12*((1/(l-nu))**13-(1/(l-nu))**7)
    fy                  = fy_bottom_boundary+fy_top_boundary
    return fx,fy


def main():
    '''======================== Random Positions ==========================='''
    '''
    for i in range(0,len(f.x_dum)):
        for j in range(i,len(f.x_dum)):
            while (np.sqrt((f.x_dum[i]-f.x_dum[j])**2+(f.y_dum[i]-f.y_dum[j])**2)<1) or (np.sqrt((f.x_dum[i]-f.l)**2+(f.y_dum[i]-f.l)**2)<1) or (np.sqrt(f.x_dum[i]**2+f.y_dum[i]**2)<1):
                print("too close together i {} j {}".format(i,j))
                f.x_dum[i]=random.random()*(f.l-1)
                f.y_dum[i]=random.random()*(f.l-1)
                if(i==j):
                    break
    '''
    for i in range(0,len(f.x_dum)):
        for j in range(i,len(f.x_dum)):
            while (np.sqrt((f.x_dum[i]-f.x_dum[j])**2+(f.y_dum[i]-f.y_dum[j])**2)<1) or (np.sqrt((f.x_dum[i]-f.l)**2+(f.y_dum[i]-f.l)**2)<1) or (np.sqrt(f.x_dum[i]**2+f.y_dum[i]**2)<1) or f.x_dum[i]<1 or f.y_dum[i]<1:
                #print("too close together i {} j {}".format(i,j))
                f.x_dum[i]=random.random()*(f.l-1)
                f.y_dum[i]=random.random()*(f.l-1)
                if(i==j):
                    break

    f.x_dum[0]=1
    f.y_dum[0]=1
    print("randomization done --- %s seconds --- changing" % (time.time()-start_time))
    #random_positions()

    r=np.zeros((f.natom,f.natom))
    u_lj_mat=np.zeros((f.natom,f.natom))

    '''======================== Verlet + Energy ==========================='''
    for n_dt in range(0,f.ndt):
        #check that everytime we restart we start from new position and origin initial cond. -cp
        f.dt_dum=array('f',[0.01*(n_dt+1)]*f.natom)
        f.ke_avg_list=array('f',[0]*f.natom)
        f.u_lj_avg_list=array('f',[0]*f.natom)
        f.u_wall_avg_list==array('f',[0]*f.natom)
        
        print("ndt: %s --- %s seconds --- changing" % (f.dt_dum[0],time.time()-start_time))
        
        for n_time in range(0,f.nstep):
            u_wall_list=array('f',[0]*f.natom)
            ke_list=array('f',[0]*f.natom)
            u_lj_list=array('f',[0]*f.natom)
            for i in range(0,f.natom):
                ''' Lennard-Jones '''
                for j in range(i,f.natom):
                    r[i][j]=((f.x_dum[i]-f.x_dum[j])**2+(f.y_dum[i]-f.y_dum[j])**2)**1/2.
                    r[j][i]=r[i][j]
                    if i==j:
                        r[i][j]=0
                        u_lj_mat[i][j]=0
                    else:
                        u_lj_mat[i][j]=lennard_jones_potential(r[i][j])
                        u_lj_mat[j][i]=u_lj_mat[i][j]
            
                ''' Each Atom's Instantaneous Energy '''
                f.u_wall[0]=np.sum(u_lj_mat[i])
                f.ke[0]=0.5*f.m_dum[i]*(f.vx_dum[i]**2+f.vy_dum[i]**2)
                f.u_lj[0]=wall_potential(f.x_dum[i],f.y_dum[i])
                
                ''' Add Energy of each Atom to List '''
                #print("u_wall: {}\t ke: {}\t u_lj: {}\t n: {} dt: {}".format(u_wall[0],ke[0],u_lj[0],n_time,dt_dum[0]))
                u_lj_list[i]=np.sum(u_lj_mat[i])
                ke_list[i]=0.5*f.m_dum[i]*(f.vx_dum[i]**2+f.vy_dum[i]**2)
                u_wall_list[i]=wall_potential(f.x_dum[i],f.y_dum[i])
                
                ''' Verlet-Method '''
                f.t_dum[i]=f.dt_dum[i]*n_time
                f.ux_dum[i]=f.vx_dum[i]+f.dt_dum[i]*f.fx_dum[i]/(2*f.m_dum[i])
                f.uy_dum[i]=f.vy_dum[i]+f.dt_dum[i]*f.fy_dum[i]/(2*f.m_dum[i])
                f.x_dum[i]=f.x_dum[i]+f.ux_dum[i]*f.dt_dum[i]
                f.y_dum[i]=f.y_dum[i]+f.uy_dum[i]*f.dt_dum[i]
                f.fx_dum[i],f.fy_dum[i]=wall_force(f.x_dum[i],f.y_dum[i])
                f.vx_dum[i]=f.ux_dum[i]+f.dt_dum[i]*f.fx_dum[i]/(2*f.m_dum[i])
                f.vy_dum[i]=f.uy_dum[i]+f.dt_dum[i]*f.fy_dum[i]/(2*f.m_dum[i])

                #triggering condition
                if n_time<1000:
                    f.t[0]=f.t_dum[i]
                    f.ux[0]=f.ux_dum[i]
                    f.uy[0]=f.uy_dum[i]
                    f.x[0]=f.x_dum[i]
                    f.y[0]=f.y_dum[i]
                    f.fx[0]=f.fx_dum[i]
                    f.fy[0]=f.fy_dum[i]
                    f.vx[0]=f.vx_dum[i]
                    f.vy[0]=f.vy_dum[i]
                    f.n[0]=n_time
                    
                    f.id[0]=i
                    f.dt[0]=f.dt_dum[i]
                    f.m[0]=f.m_dum[i]
                    
                    f.tree.Fill()
        
            if n_time<1000:
                #print("time: {}".format(t_dum[p]))
                f.u_lj_avg[0]=np.sum(u_lj_list)/(f.natom)
                f.ke_avg[0]=np.sum(ke_list)/(f.natom)
                f.u_wall_avg[0]=np.sum(u_wall_list)/(f.natom)
                f.tree2.Fill()

                f.u_lj_avg_list.append(f.u_lj_avg[0])       #not stored
                f.ke_avg_list.append(f.ke_avg[0])           #not stored
                f.u_wall_avg_list.append(f.u_wall_avg[0])   #not stored

        ''' dt change '''
        f.U_lj_time_avg     = np.sum(f.u_lj_avg_list)/len(f.u_lj_avg_list)
        f.U_wall_time_avg   = np.sum(f.u_wall_avg_list)/len(f.u_wall_avg_list)
        f.E_p_n[0]          = (f.U_lj_time_avg+f.U_wall_time_avg)/len(f.u_lj_avg_list)
        f.KE_p_n[0]         = np.sum(f.ke_avg_list)/len(f.ke_avg_list)
        f.dt[0]             = f.dt_dum[0]
        f.tree3.Fill()

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
