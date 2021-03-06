"""
*  author: carlos p
*  purpose: solve n-particle verlet
*  usage: python write_structured_md.py <natoms> <nsteps>
*
"""

import sys
import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
import numpy as np
import random
import os
from datetime import datetime
from array import array
import time

DBG=0
DBG2=0

def lj_force(xij,yij):
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2Body Separation Distance
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
    
    return fx, fy, r

def lj_potential(xij,yij):
    #lennard jones scalar potential
    #with spline compensation
    ri = 1.1080
    rc = 1.5475
    r  = np.sqrt(xij**2 + yij**2) #2-body separation distance
    if r<ri:
        if r==0:
            print("lj_potential()->avoid divide by zero")
            u = 10000
        else:
            u = (1/r)**12-2*(1/r)**6
    if ri<=r<rc:
        u = -6.129*(r-rc)**2 + -4.655*(r-rc)**3
    elif r>=rc:
        u = 0
    
    return u

class Node(object):
    def __init__(self,x=0,y=0,vx=0,vy=0,id=0,nextnode=None):
        self.id=id
        self.x=x
        self.y=y
        self.cellx=0
        self.celly=0
        self.fx=[]
        self.fy=[]
        self.p=[]
        self.ulj=[]
        self.r=[]
        self.vx=vx
        self.vy=vy
        self.ux=0
        self.uy=0
        self.m=0
        self.dt=0
        self.iter_issue=None
        self.next_node=None
        self.prev_node=None
    
    def set_next_node(self,n):
        self.next_node=n
    
    def get_next_node(self):
        return self.next_node

class LinkedList(object):
    def __init__(self,r=None): #,r=None
        self.root=None
        self.natoms=0
        self.alpha=0
        self.temp=0
        self.dt=0
        self.kb=1
        self.a=0
        self.b=0
        self.c=0
        self.d=0
        self.L=0
        self.nelements=0
        self.nx=0 #x cell element no.
        self.ny=0 #y cell element no.
    
    def add(self,node): #reexamine this routine nextnode setting
        #print("\tAdd node id={}".format(node.id))
        if node!=None:
            if self.root is None:
                node.next_node=None #self.root #first
                node.prev_node=None
                self.root=node
                #print("LinkedList()->add() first x={} yx={}".format(node.x,node.y))
                
            elif self.root is not None and self.root.next_node is None: #second
                first_node=self.root
                first_node.next_node=None
                first_node.prev_node=node #old node point back
                node.next_node=first_node #new node point forward
                node.prev_node=None #new node point to nothing
                self.root=node #start at inserted node
                #print("LinkedList()->add() secondx={} yx={}".format(node.x,node.y))
                
            elif self.root is not None and self.root.next_node is not None:
                first_node=self.root
                first_node.prev_node=node
                node.next_node=first_node
                node.prev_node=None
                self.root=node
        
            if node==node.next_node:
                print("We FUCKed up id={} id2={}".format(node.id,node.next_node.id))
                #return False
                exit()
                return False
            else:
                return True

    def removeNode(self,node):
        #return node to be added to other cell of grid
        
        #linked_node=self.root
        
        #remove node from sequence and return the address
        #while linked_node is not None:
        #if linked_node.id is node.id: #select node
        if node!=None:
            #case0: linked node is only one in cell
            if node.next_node is None and node.prev_node is None:
                self.root=None
                if DBG: print("IM STUCK 1")
                return node
            
            #case1: linked node is in middle of list
            if node.next_node is not None and node.prev_node is not None:
                first_node=node.next_node
                second_node=node.prev_node
                
                first_node.prev_node=second_node #point back at prev node
                second_node.next_node=first_node
                if DBG: print("IM STUCK 2")
                #return linked_node
                return node
            
            #case2: linked node is root
            if node.prev_node is None and node.next_node is not None:
                first_node=node.next_node
                first_node.prev_node=None
                self.root=first_node
                if DBG: print("IM STUCK 3")
                #return linked_node
                return node
            
            #case3: linked node is first
            if node.next_node is None and node.prev_node is not None:
                second_node=node.prev_node
                second_node.next_node=None
                if DBG: print("IM STUCK 4")
                #return linked_node
                return node
    
        if DBG: print("removeNode->END")
        #linked_node=linked_node.next_node

    def getSize(self):
        return self.natoms

    def getEnergies(self,tn):
        node=self.root
        ekv=0
        eku=0
        ulj=0
        counter=0
        while node is not None:
            ekv+=0.5*node.m*(node.vx**2+node.vy**2)
            eku+=0.5*node.m*(node.ux**2+node.uy**2)
            ulj+=np.sum(node.ulj)
            #print("id={} vx={} ekv={} eku={} ulj={}".format(node.id,node.vx,ekv,eku,ulj))
            if tn==0 and DBG2:
                print("tn={} vx={} vy={} ux={} uy={}".format(tn,node.vx,node.vy,node.ux,node.uy))
            
            node=node.next_node
        
        return ekv,eku,ulj

    def updateInteractionForces(self,linked_list_b,n):
        node_a=self.root #root of linked list
        node_b=linked_list_b.root
        list_b=linked_list_b
        label="Collision: (id={})=[{},{}] (id={})=[{},{}] r={} fx={:.2f} fy={:.2f}"
        same_node=False
        countera=0

	#return values
	error=False
	p=0
        while node_a!=None:
            if countera>50 or error==True:
                label="Linked()..Error: ForceInter a idx={} nodes[prev-{} cur-{} next-{} nextX={}"
                print(label.format(n,node_a.prev_node,node_a,node_a.next_node,node_a.next_node.x))

                error=True
		break

	    else:
		error=False
            
            countera+=1
            counterb=0
            while node_b!=None:
                if counterb>50:
                    label="Linked()..Error: ForceInter b idx={} nodes[prev={} cur={} next={}->x={}"
                    print(label.format(n,node_b.prev_node,node_b,node_b.next_node,node_b.next_node.x))
                    error=True
		    break

		else:
		    error=False
                
                if node_a.id==node_b.id:
                    if DBG: print("ForceInt: same node collision [{},{}]".format(node_a.cellx,node_b.celly))
                    node_b=node_b.next_node
                    break
                
                xij=node_a.x-node_b.x
                yij=node_a.y-node_b.y
                if xij>0.5*self.L:   xij-=self.L
                if xij<=-0.5*self.L: xij+=self.L
                if yij>0.5*self.L:   yij-=self.L
                if yij<=-0.5*self.L: yij+=self.L
                xji=-xij
                yji=-yij

                fxij,fyij,r=lj_force(xij,yij)
                fxji=-fxij
                fyji=-fyij

                if abs(fxij)>10000 or (abs(fyij)>10000):
                    node_a.iter_issue=n
                    label="High Force at time={} id_a={} id_b={} - xij={}"
                    if DBG: print(label.format(n,node_a.id,node_b.id,xij))
                    break #end time iter'''
                
                ulj=lj_potential(xij,yij)
                if r<1.54 and (fxij!=0 or fyij!=0) and 1:
                    if DBG2: print("====================================================")
                    if DBG2: print(label.format(node_a.id,self.nx,self.ny
                                       ,node_b.id,list_b.nx,list_b.ny,r,fx,fy))
                
                node_a.fx.append(fxij)
                node_a.fy.append(fyij)
                node_a.ulj.append(ulj)
                node_b.fx.append(fxji)
                node_b.fy.append(fyji)
                node_b.ulj.append(ulj)
                node_a.r.append(r)

                #pressure algorithm
                pij=0.5*(xij*fxij + yij*fyij)
                pji=0.5*(xji*fxji + yji*fyji)
		p+=pij                

                node_a.p.append(pij)
                node_b.p.append(pji)

		#next iteration setup
		counterb+=1
                node_b=node_b.next_node
                if node_b is None:
                    break

	    #next iteration setup
            node_a=node_a.next_node

	return error,p

    def clearNodeForces(self):
        node_a=self.root #linked list root
        counter=0
        while node_a is not None:
            #print("clear force id: {} - force: {}".format(node_a.id,node_a.fx))
            del node_a.fx[:]
            del node_a.fy[:]
            del node_a.p[:] #clear pressure contribution also
            del node_a.ulj[:] #clear potentials also
            node_a=node_a.next_node
            counter+=1
            if counter>50:
                print("clearNodeForces")
                return False

        if DBG: print("Linked->clearNode->END")

    def isEmpty(self):
        if self.root is None:
            return None

    def updateListedPositions(self,tn):
        title="LinkedList()->updateListedPositions particle id={}"
        node=self.root
        
        if node is not None:
            i2m=1/(2.*node.m)
            self.a=(1-self.alpha*self.dt*i2m )/(1+self.alpha*self.dt*i2m )
            self.b=1/(1+self.alpha*self.dt*i2m )
            self.c=self.a/np.sqrt(self.b)
            self.d=np.sqrt(self.b)
            self.updateflag=False

        # Langevin Motion Algorithm
        while node is not None:
            if tn==0 and DBG2:
                print("uPos tn={} vx={} vy={} ux={} uy={}".format(tn,node.vx,node.vy,node.ux,node.uy))
            
            fx=np.sum(node.fx)
            fy=np.sum(node.fy)
            
            xi1=random.uniform(0,1)
            if xi1==0: xi1=1-xi1
                        
            xi2=random.uniform(0,1)
            if xi2==0: xi2=1-xi2
            
            eta1=np.sqrt(-2*np.log(xi1))*np.cos(2*np.pi*xi2)
            eta2=np.sqrt(-2*np.log(xi1))*np.sin(2*np.pi*xi2)
            beta1=np.sqrt(2*self.alpha*self.kb*self.temp*self.dt)*eta1
            beta2=np.sqrt(2*self.alpha*self.kb*self.temp*self.dt)*eta2
            # turfv
            node.t=self.dt*(tn+1)
            node.ux=self.d*(node.vx+(self.dt*fx+beta1)*i2m )
            node.uy=self.d*(node.vy+(self.dt*fy+beta2)*i2m )
            node.x=node.x+self.d*self.dt*node.ux
            node.y=node.y+self.d*self.dt*node.uy
            
            if node.id is 0 and DBG2:
                test="a={} b={} c={} d={} i2m={} dt={} beta={} x={}"
                print(test.format(self.a,self.b,self.c,self.d,i2m,self.dt,beta1,node.x))
            
            node.vx=self.c*node.ux+(self.dt*fx+beta1)*i2m
            node.vy=self.c*node.uy+(self.dt*fy+beta2)*i2m
            
            #reenter box boundary
            if node.x>=self.L: node.x-=self.L
            if node.x<0: node.x+=self.L
            if node.y>=self.L: node.y-=self.L
            if node.y<0: node.y+=self.L
            
            node=node.next_node

    #def save2root(self):
    def reset(self):
        self.root=None
        self.natoms=0
        self.alpha=0
        self.temp=0
        self.dt=0
        self.kb=1
        self.a=0
        self.b=0
        self.c=0
        self.d=0
        self.L=0
        self.nelements=0
        self.nx=0 #x cell element no.
        self.ny=0 #y cell element no.
    
    
    def printList(self):
        node = self.root
        while node is not None:
            listStr="LinkedList->printList() (id={})=[{},{}]=[{:.2f},{:.2f}] vx={} vy={} dt={}"
            print(listStr.format(node.id,node.cellx,node.celly,node.x,node.y,node.vx,node.vy,node.dt))
            node = node.next_node

class LinkedListMat(object):
    def __init__(self):
        self.pixel_root=None
        self.natoms=0
        self.alpha=0
        self.temp=0
        self.ndt=0
        self.vx0=0
        self.vy0=0
        self.nx=0       #nx elements
        self.ny=0       #ny elements
        self.L=20
    
    def addGrid(self,nx,ny):
        linked_list=LinkedList()
        self.pixel_root=[[linked_list for j in range(nx)] for i in range(ny)]
        self.nx=nx
        self.ny=ny
        self.natoms=0
        #grid of linked list
        for i in range(nx):
            for j in range(ny):
                linked_list=LinkedList()
                linked_list.nx=i
                linked_list.ny=j
                self.pixel_root[i][j]=linked_list
    
    def addAtoms(self,natoms):
        #add n atoms
        linked_list = LinkedList()
        xj=[]
        yj=[]
        
        for nn in range(0,natoms):
            node=Node()
            node.x=random.random()*(self.L)
            node.y=random.random()*(self.L)
            for nnn in range(len(xj)):
                xij=node.x-xj[nnn]
                yij=node.y-yj[nnn]
                hyp=(xij**2 + yij**2)**0.5
                counter=0
                while hyp<2 and nn is not nnn:
                    if DBG2: print("adjuments to id={} and previous id={}".format(nn,nnn))
                    node.x=random.random()*(self.L)
                    node.y=random.random()*(self.L)
                    xij=node.x-xj[nnn]
                    yij=node.y-yj[nnn]
                    hyp=(xij**2 + yij**2)**0.5
                    if counter>100:
                        print("addAtoms:Error too many tries")
                        exit(1)
                    
                    counter+=1
            
            node.cellx=int(node.x*(self.nx/(1.*self.L) ))
            node.celly=int(node.y*(self.ny/(1.*self.L) ))
            node.vx=2
            node.vy=3
            node.ux=2
            node.uy=3
            node.m=1
            node.id=nn
            #print("id:{}".format(node.id))
            self.pixel_root[node.cellx][node.celly].add(node) #should this be =?
            #return root to pixel root
            self.natoms+=1
            xj.append(node.x)
            yj.append(node.y)
            
        self.vx0=node.vx
        self.vy0=node.vy

    def setParameters(self,alpha,temp,dt):
        #set grid attributes of linked list
        for i in range(self.nx):
            for j in range(self.ny):
                ll=self.pixel_root[i][j]
                if ll is not None:
                    ll.alpha=alpha
                    ll.temp=temp
                    ll.dt=dt
                    ll.L=self.L

    def updateForces(self,n):
        #update forces within grid
        nx=self.nx-1 #-1 correction
        ny=self.ny-1 #-1 correction
        p=0
        for n_x in range(0,nx+1):
            for n_y in range(0,ny+1):
                cell=self.pixel_root[n_x][n_y]
                node=cell.root
                cell2=None
                cell3=None
                cell4=None
                cell5=None
                #determine neighboring cells
                if node is not None:
                    #top condition DONE
                    if n_x<nx  and n_x>=1 and n_y==ny:
                        c2=self.pixel_root[n_x+1][ny]
                        c3=self.pixel_root[n_x+1][0]
                        c4=self.pixel_root[n_x][0]
                        c5=self.pixel_root[n_x-1][0]
                        
                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                    
                    #right corner condition DONE
                    if n_x==nx and n_y==ny:
                        c2=self.pixel_root[0][ny]        #right
                        c3=self.pixel_root[0][0]         #top right
                        c4=self.pixel_root[nx][0]        #top
                        c5=self.pixel_root[nx-1][0]      #top left
                    
                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                    
                    #left corner condition DONE
                    if n_x==0  and n_y==ny:
                        c2=self.pixel_root[n_x+1][ny]    #right
                        c3=self.pixel_root[n_x+1][0]     #top right
                        c4=self.pixel_root[n_x][0]       #top
                        c5=self.pixel_root[nx][0]        #top left
                    
                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                        
                    #right side condition DONE
                    if n_x==nx and n_y<ny:
                        c2=self.pixel_root[0][n_y]
                        c3=self.pixel_root[0][n_y+1] #set ny=0 - out of bounds
                        c4=self.pixel_root[n_x][n_y+1] #set ny=0 - out of bounds
                        c5=self.pixel_root[n_x-1][n_y+1] #set ny=0 - out of bounds
                    
                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                    
                    #left side condition
                    if n_x==0  and n_y<ny:
                        c2=self.pixel_root[n_x+1][n_y]
                        c3=self.pixel_root[n_x+1][n_y+1]
                        c4=self.pixel_root[n_x][n_y+1]
                        c5=self.pixel_root[nx][n_y+1]

                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                            
                    #normal condition
                    if n_x<nx and n_x>=1 and n_y<ny and n_y>=0:
                        c2=self.pixel_root[n_x+1][n_y]
                        c3=self.pixel_root[n_x+1][n_y+1]
                        c4=self.pixel_root[n_x][n_y+1]
                        c5=self.pixel_root[n_x-1][n_y+1]

                        if c2 is not None: cell2=c2     #right
                        if c3 is not None: cell3=c3     #right top
                        if c4 is not None: cell4=c4     #top
                        if c5 is not None: cell5=c5     #top left
                    
                    node = node.next_node
                    
                    done1=True
                    done2=True
                    done3=True
                    done4=True
                    done5=True
                    #calculate force
                    if cell is not None: err1,p1=cell.updateInteractionForces(cell,n)
                    if cell2 is not None: err2,p2=cell.updateInteractionForces(cell2,n)
                    if cell3 is not None: err3,p3=cell.updateInteractionForces(cell3,n)
                    if cell4 is not None: err4,p4=cell.updateInteractionForces(cell4,n)
                    if cell5 is not None: err5,p5=cell.updateInteractionForces(cell5,n)

                    if err1 is True: print("cell->cell interaction error")
                    if err2 is True: print("cell->cell2 interaction error")
                    if err3 is True: print("cell->cell3 interaction error")
                    if err4 is True: print("cell->cell4 interaction error")
                    if err5 is True: print("cell->cell5 interaction error")

                    p=p+p1+p2+p3+p4+p5

        p=p/self.natoms
        if DBG: print("IntForce p={} natoms={}".format(p,self.natoms))
        return p

    def clearForces(self):
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                node=cell.root
                if node is not None:
                    done=cell.clearNodeForces()
                    if done is False:
                        print("ClearForces-natoms={}".format(self.natoms))
                        exit()

    def reset(self):
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                #cell=None
                node=cell.root
                
                self.natoms=0
                counter=0
                while node!=None:
                    node=cell.removeNode(node)
                    counter+=1
                    if node!=None:
                        node=node.next_node
                    else:
                        break
                    
                    if counter>50:
                        print("Mat->reset next={}".format(node.next))
                        exit()

                if DBG: print("Mat->Reset()->END")

    def updatePositions(self,tn):
        kev_mat=[]
        keu_mat=[]
        u_lj_mat=[]
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                if cell is not None:
                    cell.updateListedPositions(tn)

    def getAvgEnergy(self,tn):
        kev_mat=[]
        keu_mat=[]
        ulj_mat=[]
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                if cell is not None:
                    #get energies
                    kev,keu,ulj=cell.getEnergies(tn)
                    kev_mat.append(kev)
                    keu_mat.append(keu)
                    ulj_mat.append(ulj)

        ekv=np.sum(kev_mat)/self.natoms
        eku=np.sum(keu_mat)/self.natoms
        epn=np.sum(ulj_mat)/self.natoms
        if DBG: print("linkMat->getAvgEnergy(): natoms={}".format(self.natoms))
        return ekv,eku,epn

    def updateCells(self,tn):
        #print("updateCells")
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                linked_list=self.pixel_root[nx][ny]
                out_of_bounds=False
                
                if linked_list is not None:
                    node=linked_list.root
                    
                    counter=0
                    while node!=None:
                        cur_nx=node.cellx
                        cur_ny=node.celly
                        new_nx=int(node.x*(self.nx/(1.*self.L) ))
                        new_ny=int(node.y*(self.ny/(1.*self.L) ))
                        
                        if cur_nx!=new_nx or cur_ny!=new_ny:
                            if node.x<self.L and node.x>0 and node.y<self.L and node.y>0:
                                label="Update id:{} Cell Assignement [{},{}]->[{},{}]"
                                if DBG: print(label.format(node.id,cur_nx,cur_ny,new_nx,new_ny))
                                node.cellx=new_nx
                                node.celly=new_ny
                                
                                node=linked_list.removeNode(node) #node address
                                self.pixel_root[new_nx][new_ny].add(node)
                                node=node.next_node #move to next node
                            
                            else:
                                counter+=1
                                node.x=random.random()*(self.L-1)
                                node.y=random.random()*(self.L-1)
                                
                                if counter>50:
                                    label="im stuck fx: {} fy: {} time-iter={} nx={} ny={}"
                                    print(label.format(np.sum(node.fx),
                                                       np.sum(node.fy),
                                                       tn,
                                                       new_nx,
                                                       new_ny) )
                                    exit()
                    
                        else:
                            if DBG: print("updateCells Ctr={} NextNode={}".format(counter,node.next_node))
                            counter+=1
                            node=node.next_node
                            if counter>50 or node is None:
                                break
                                #exit()

                        if DBG: print("updateCells->END")

    def printList(self):
        print("PrintList====================================")
        for n_x in range(0,self.nx):
            for n_y in range(0,self.ny):
                linked_list = self.pixel_root[n_x][n_y]
                #linked_list.printList()
                if linked_list is not None:
                    node=linked_list.root
                    if node is not None:
                        listStr="LinkedList->printList() (id={})=[{},{}]=[{:.2f},{:.2f}] vx={} vy={} dt={}"
                        print(listStr.format(node.id,
                                             node.cellx,
                                             node.celly,
                                             node.x,
                                             node.y,
                                             node.vx,
                                             node.vy,
                                             node.dt))

filename="data/celldat_{}natom_{}nsteps_vx{:.0f}_vy{:.0f}_{}ndt.root"
def main():
    #control variables - eventually feed through commandline args -cp
    nrep=1
    natomstep=1
    natommax=int(sys.argv[1]) #30
    nx=10
    ny=10
    nsteps=int(sys.argv[2])
    alphaList=[1]
    tempList=[10]
    dtstep=0.001
    dtmax =0.001
    ndt=int(dtmax/dtstep)
    
    #general variables
    file=TFile(filename.format(natommax,nsteps,2,3,ndt),"recreate")
    nreps=array('i',[nrep])
    pressure=array('f',[0]*nreps[0])
    alpha_id=array('i',[0])
    temp_id=array('i',[0])
    dt_id=array('i',[0])
    alpha=array('f',[0])
    temp=array('f',[0])
    dt=array('f',[0])
    #study variables
    natom=array('i',[0]*nreps[0])
    runtime=array('f',[0]*nreps[0])
    ekvAvg=array('f',[0]*nreps[0])
    ekuAvg=array('f',[0]*nreps[0])
    epvAvg=array('f',[0]*nreps[0])
    epuAvg=array('f',[0]*nreps[0])
    
    tree=TTree("averageEnergy", "energy analysis")
    tree.Branch('nreps',nreps,'nreps/I')
    tree.Branch('pressure',pressure,'pressure[nreps]/F')
    tree.Branch('alpha_id',alpha_id,'alpha_id/I')
    tree.Branch('alpha',alpha,'alpha/F')
    tree.Branch('temp_id',temp_id,'temp_id/I')
    tree.Branch('temp',temp,'temp/F')
    tree.Branch('dt_id',dt_id,'dt_id/I')
    tree.Branch('dt',dt,'dt/F')
    tree.Branch('ekvAvg',ekvAvg,'ekvAvg[nreps]/F')
    tree.Branch('ekuAvg',ekuAvg,'ekuAvg[nreps]/F')
    tree.Branch('epvAvg',epvAvg,'epvAvg[nreps]/F')
    tree.Branch('epuAvg',epuAvg,'epuAvg[nreps]/F')
    tree.Branch('runtime',runtime,'runtime[nreps]/F')
    tree.Branch('natom',natom,'natom[nreps]/I')

    #generate grid
    linkedListMatrix=LinkedListMat()
    linkedListMatrix.addGrid(nx,ny)
    
    for n_alpha in range(0,len(alphaList)):
        alpha[0]=alphaList[n_alpha]
        
        for n_temp in range(0,len(tempList)):
            temp[0]=tempList[n_temp]
            
            for n_atom in range(0,natommax):
                start_time=time.time()

                for n_dt in range(0,1):
                    dt[0]=(1)*dtstep #n_dt+

                    for n_rep in range(0,nreps[0]):
                        ekv_sum=0
                        eku_sum=0
                        ulj_sum=0
                        p_sum=0
                        natom[n_rep]=(n_atom+1)*natomstep

                        linkedListMatrix.reset()
                        linkedListMatrix.addAtoms(natom[n_rep])
                        linkedListMatrix.setParameters(alpha[0],temp[0],dt[0])

                        for n_time in range(0,nsteps):
                            p_sum+=linkedListMatrix.updateForces(n_time)
                            linkedListMatrix.updatePositions(n_time)
                            ekv,eku,ulj=linkedListMatrix.getAvgEnergy(n_time)
                            ekv_sum+=ekv
                            eku_sum+=eku
                            ulj_sum+=ulj
                            linkedListMatrix.clearForces()
                            linkedListMatrix.updateCells(n_time)

                        #average system energy
                        ekvAvgTest=ekv_sum/nsteps
                        ekuAvgTest=eku_sum/nsteps
                        uljAvgTest=ulj_sum/nsteps
                        pavgTest=p_sum/nsteps

                        #save results2root
                        pressure[n_rep]=pavgTest
                        ekvAvg[n_rep]=ekvAvgTest
                        ekuAvg[n_rep]=ekuAvgTest
                        epvAvg[n_rep]=ekvAvgTest+uljAvgTest
                        epuAvg[n_rep]=ekuAvgTest+uljAvgTest
                        
                        #id info
                        alpha_id[0]=0
                        temp_id[0]=0
                        dt_id[0]=n_dt

                        runtime[0]=float(time.time()-start_time)
                        tree.Fill()

                        #print results
                        timelabel1="time={:.2f} natom={} "
                        label1="dt={:.3f}: p={:.1f} ekv={:.2f} eku={:.2f} ulj={:.2f} "
                        label2="[alpha={:.1f} temp={:.1f}]"
                        print((timelabel1+label1+label2).format(runtime[0],
                                                                natom[n_rep],
                                                                dt[0],
                                                                pressure[n_rep],
                                                                ekvAvg[n_rep],
                                                                ekuAvg[n_rep],
                                                                epuAvg[n_rep],
                                                                alpha[0],
                                                                temp[0]))
                                                                
                        #linkedListMatrix.printList()
                        
                    #write
                    tree.Write()

    file.Close()
                                     

if __name__ == "__main__":
    main()


