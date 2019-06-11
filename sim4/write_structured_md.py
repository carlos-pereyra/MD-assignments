"""
*  author: carlos p
*  purpose: solve n-particle verlet
*
"""

import ROOT
from ROOT import TTree,TFile,gROOT
from ROOT import TCanvas,TGraph,TPad,TBrowser
import numpy as np
import random
import os
from datetime import datetime
from array import array

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
        self.ulj=[]
        self.r=[]
        self.vx=vx
        self.vy=vy
        self.m=0
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
    
    def add(self,node):
        if node is not None:
            if self.root is None:
                node.next_node=self.root #first
                node.prev_node=None
                self.root=node
                #print("LinkedList()->add() first x={} yx={}".format(node.x,node.y))
            elif self.root is not None and self.root.prev_node is None: #second
                first_node=self.root
                first_node.prev_node=node #old node point back
                node.next_node=first_node #new node point forward
                node.prev_node=None #new node point to nothing
                self.root=node #start at inserted node
                #print("LinkedList()->add() secondx={} yx={}".format(node.x,node.y))
            elif self.root is not None and self.root.prev_node is not None:
                first_node=self.root
                first_node.prev_node=node
                node.next_node=first_node
                node.prev_node=None
                self.root=node

    def removeNode(self,node):
        linked_node=self.root
        
        #remove node from sequence and return the address
        while linked_node is not None:
            if linked_node.id is node.id: #select node
                #case0: linked node is only one in cell
                if linked_node.next_node is None and linked_node.prev_node is None:
                    self.root=None
                    if DBG: print("IM STUCK 1")
                    return linked_node
                
                #case1: linked node is in middle of list
                if linked_node.next_node is not None and linked_node.prev_node is not None:
                    first_node=linked_node.next_node
                    second_node=linked_node.prev_node
                    
                    first_node.prev_node=second_node #point back at prev node
                    second_node.next_node=first_node
                    if DBG: print("IM STUCK 2")
                    return linked_node
                
                #case2: linked node is root
                if linked_node.prev_node is None and linked_node.next_node is not None:
                    first_node=linked_node.next_node
                    first_node.prev_node=None
                    self.root=first_node
                    if DBG: print("IM STUCK 3")
                    return linked_node
                
                #case3: linked node is first
                if linked_node.next_node is None and linked_node.prev_node is not None:
                    second_node=linked_node.prev_node
                    second_node.next_node=None
                    if DBG: print("IM STUCK 4")
                    return linked_node
                
            linked_node=linked_node.next_node

    def getSize(self):
        return self.natoms

    def getEnergies(self):
        node=self.root
        ekv=0
        eku=0
        ulj=0
        while node is not None:
            ekv+=0.5*(node.vx*node.vx+node.vy*node.vy)
            eku+=0.5*(node.ux*node.ux+node.uy*node.uy)
            ulj+=np.sum(node.ulj)
            #print("id={} vx={} ekv={} eku={} ulj={}".format(node.id,node.vx,ekv,eku,ulj))
            node=node.next_node
        
        return ekv,eku,ulj

    def updateInteractionForces(self,linked_list_b,n):
        node_a=self.root #root of linked list
        node_b=linked_list_b.root
        list_b=linked_list_b
        label="Collision: (id={})=[{},{}] (id={})=[{},{}] r={} fx={:.2f} fy={:.2f}"
        same_node=False
        #empty stored node forces
        while node_a is not None:
            while node_b is not None:
                if node_a.id is node_b.id and node_a is not 0:
                    #print("same cell collision [{},{}]".format(node_a.cellx,node_b.celly))
                    same_node=True
                    node_b=node_b.next_node
                    continue
                
                xij=node_a.x-node_b.x
                yij=node_a.y-node_b.y
                
                if xij>0.5*self.L:   xij-=self.L
                if xij<=-0.5*self.L: xij+=self.L
                if yij>0.5*self.L:   yij-=self.L
                if yij<=-0.5*self.L: yij+=self.L
                
                fx,fy,r=lj_force(xij,yij)
                if abs(fx)>2000 or (abs(fx)>2000):
                    node_a.iter_issue=n
                    print("High Force at time={} id_a={} id_b={} - xij={}".format(n,node_a.id,node_b.id,xij))
                    break #end time iter
                
                ulj=lj_potential(xij,yij)
                if r<1.54 and (fx!=0 or fy!=0) and 1:
                    if DBG2: print("====================================================")
                    if DBG2: print(label.format(node_a.id,self.nx,self.ny
                                       ,node_b.id,list_b.nx,list_b.ny,r,fx,fy))
                
                node_a.fx.append(fx)
                node_a.fy.append(fy)
                node_a.ulj.append(ulj)
                node_b.fx.append(-fx)
                node_b.fy.append(-fy)
                node_b.ulj.append(ulj)
                node_a.r.append(r)
                node_b=node_b.next_node
            
            if same_node is True:
                break

            node_a=node_a.next_node

    def clearNodeForces(self):
        node_a=self.root #linked list root
        while node_a is not None:
            #print("clear force id: {} - force: {}".format(node_a.id,node_a.fx))
            del node_a.fx[:]
            del node_a.fy[:]
            del node_a.ulj[:] #clear potentials also
            node_a=node_a.next_node
    

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
            #print(title.format(node.id))
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
            node.t=node.dt*(tn+1)
            node.ux=self.d*(node.vx+(node.dt*fx+beta1)*i2m )
            node.uy=self.d*(node.vy+(node.dt*fy+beta2)*i2m )
            node.x=node.x+self.d*node.dt*node.ux
            node.y=node.y+self.d*node.dt*node.uy
            node.vx=self.c*node.ux+(self.dt*fx+beta1)*i2m
            node.vy=self.c*node.uy+(self.dt*fy+beta2)*i2m
            
            #reenter box boundary
            if node.x>=self.L: node.x-=self.L
            if node.x<0: node.x+=self.L
            if node.y>=self.L: node.y-=self.L
            if node.y<0: node.y+=self.L
            
            node=node.next_node

    #def save2root(self):

    def printList(self):
        node = self.root
        while node is not None:
            listStr="LinkedList->printList() (id={})=[{},{}]=[{:.2f},{:.2f}]"
            print(listStr.format(node.id,node.cellx,node.celly,node.x,node.y))
            node = node.next_node

class LinkedListMat(object):
    def __init__(self):
        self.pixel_root=None
        self.nx=0       #nx elements
        self.ny=0       #ny elements
        self.L=20
        self.natoms=0
    
    def addGrid(self,nx,ny):
        linked_list=LinkedList()
        self.pixel_root=[[linked_list for j in range(nx)] for i in range(ny)]
        self.nx=nx
        self.ny=ny
        #grid of linked list
        for i in range(nx):
            for j in range(ny):
                linked_list=LinkedList()
                linked_list.nx=i
                linked_list.ny=j
                self.pixel_root[i][j]=linked_list
    
    def addAtoms(self,n):
        #add n atoms
        linked_list = LinkedList()
        xj=[]
        yj=[]
        
        for nn in range(n):
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
            
            node.cellx  = int(node.x*(self.nx/(1.*self.L) ))
            node.celly  = int(node.y*(self.ny/(1.*self.L) ))
            node.dt     = 0.001
            node.vx     = 2
            node.vy     = 3
            node.m      = 1
            node.id     = nn
            self.pixel_root[node.cellx][node.celly].add(node)
            self.natoms += 1
            xj.append(node.x)
            yj.append(node.y)

    def setParameters(self,alpha,temp,dtstep):
        #set grid attributes of linked list
        for i in range(self.nx):
            for j in range(self.ny):
                ll=self.pixel_root[i][j]
                if ll is not None:
                    ll.alpha=alpha
                    ll.temp=temp
                    ll.dtstep=dtstep
                    ll.L=self.L

    def updateForces(self,n):
        #update forces within grid
        nx=self.nx-1 #-1 correction
        ny=self.ny-1 #-1 correction
        for n_x in range(0,nx+1):
            for n_y in range(0,ny+1):
                cell=self.pixel_root[n_x][n_y]
                node=self.pixel_root[n_x][n_y].root
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
                    
                    #calculate force
                    if cell is not None: cell.updateInteractionForces(cell,n)
                    if cell2 is not None: cell.updateInteractionForces(cell2,n)
                    if cell3 is not None: cell.updateInteractionForces(cell3,n)
                    if cell4 is not None: cell.updateInteractionForces(cell4,n)
                    if cell5 is not None: cell.updateInteractionForces(cell5,n)

    def clearForces(self):
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                node=cell.root
                if node is not None:
                    #print("Node {}: F={}".format(node.id,node.fx))
                    #node=node.next_node
                    cell.clearNodeForces()

    def reset(self):
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                node=cell.root
                while node is not None:
                    cell.removeNode(node)
                    node=node.next_node

    def updatePositions(self,tn):
        kev_mat=[]
        keu_mat=[]
        u_lj_mat=[]
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                if cell is not None:
                    cell.updateListedPositions(tn)

    def getAvgEnergy(self):
        kev_mat=[]
        keu_mat=[]
        ulj_mat=[]
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                cell=self.pixel_root[nx][ny]
                if cell is not None:
                    #get energies
                    kev,keu,ulj=cell.getEnergies()
                    kev_mat.append(kev)
                    keu_mat.append(keu)
                    ulj_mat.append(ulj)

        ekv=np.sum(kev_mat)/self.natoms
        eku=np.sum(keu_mat)/self.natoms
        epn=np.sum(ulj_mat)/self.natoms
        return ekv,eku,epn

    def updateCells(self):
        for nx in range(0,self.nx):
            for ny in range(0,self.ny):
                linked_list=self.pixel_root[nx][ny]
                out_of_bounds=False
                if linked_list is not None:
                    node=linked_list.root
                    
                    counter=0
                    while node is not None:
                        oldcellx=node.cellx
                        oldcelly=node.celly
                        newcellx=int(node.x*(self.nx/(1.*self.L) ))
                        newcelly=int(node.y*(self.ny/(1.*self.L) ))
                        
                        '''
                        if node.x>self.L or node.y>self.L: #within bounds
                            node.x=random.random()*(self.L)
                            node.y=random.random()*(self.L)
                            newcellx=int(node.x*(self.nx/(1.*self.L)))
                            newcelly=int(node.y*(self.ny/(1.*self.L) ))
                            if 1: print("newX:{} X:{} newY:{} Y:{}".format(newcellx,node.x,newcelly,node.y))
                                #exit()
                            linked_list.removeNode(node)
                            self.pixel_root[newcellx][newcelly].add(node)'''
                        
                        counter+=1
                        if counter>50:
                            print("im stuck fx: {} fy: {} time-iter={}".format(np.sum(node.fx),np.sum(node.fy),node.iter_issue) )
                            exit()
                        
                        if node.x>self.L or node.x<0 or node.y>self.L or node.y<0:
                            node.x=random.random()*(self.L-1)
                            node.y=random.random()*(self.L-1)
                            newcellx=int(node.x*(self.nx/(1.*self.L)))
                            newcelly=int(node.y*(self.ny/(1.*self.L) ))
                            node.cellx=newcellx
                            node.celly=newcelly
                            if 1: print("try={} newX:{} X:{} newY:{} Y:{} fx={} fy={}".format(counter, newcellx, node.x, newcelly, node.y, np.sum(node.fx), np.sum(node.fx) ))
                            
                            linked_list.removeNode(node)
                            self.pixel_root[newcellx][newcelly].add(node)
                            out_of_bounds=True
                        
                        elif (out_of_bounds is not True) and (oldcellx!=newcellx or oldcelly!=newcelly):
                            node.cellx=newcellx
                            node.celly=newcelly
                            linked_list.removeNode(node)
                            label="Update id:{} Cell Assignement [{},{}]->[{},{}]"
                            if DBG: print(label.format(node.id,oldcellx,oldcelly,newcellx,newcelly))
                            self.pixel_root[newcellx][newcelly].add(node)
                            out_of_bounds=False
                        
                        if out_of_bounds is not True:
                            node=node.next_node

    def printList(self):
        print("PrintList====================================")
        for n_x in range(0,self.nx):
            for n_y in range(0,self.ny):
                linked_list = self.pixel_root[n_x][n_y]
                linked_list.printList()

def main():
    n=21
    nx=10
    ny=10
    nsteps=int(1E3)
    
    alpha=1
    
    temp=1
    
    dtstep=0.001
    dtmax =0.01
    ndt=int(dtmax/dtstep)
    
    linkedListMatrix=LinkedListMat()
    linkedListMatrix.addGrid(nx,ny)
    linkedListMatrix.printList()
    
    for n_dt in range(0,ndt):
        ekv_sum=0
        eku_sum=0
        ulj_sum=0
        dt=(n_dt+1)*dtstep
        linkedListMatrix.reset()
        linkedListMatrix.addAtoms(n) #autoadd to cells
        linkedListMatrix.setParameters(alpha,temp,dt) #comes after addGrid()
        for n_time in range(0,nsteps):
            linkedListMatrix.updateForces(n_time)
            linkedListMatrix.updatePositions(n_time)
            ekv,eku,ulj=linkedListMatrix.getAvgEnergy()
            ekv_sum+=ekv
            eku_sum+=eku
            ulj_sum+=ulj
            linkedListMatrix.updateCells()
            linkedListMatrix.clearForces()

        ekvAvg=ekv_sum/nsteps
        ekuAvg=eku_sum/nsteps
        uljAvg=ulj_sum/nsteps
        
        label="{} ndt out of {}: ekv={:.2f} eku={:.2f} ulj={:.2f} [dt={:.3f}, alpha={:.1f}, temp={:.1f}]"
        print(label.format(n_dt,ndt,ekvAvg,ekuAvg,uljAvg,dt,alpha,temp))
        linkedListMatrix.reset()

    linkedListMatrix.printList()

if __name__ == "__main__":
    main()
