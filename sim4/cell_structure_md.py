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


class Node(object):
    def __init__(self,x=0,y=0,vx=0,vy=0,id=0,nextnode=None):
        self.id=id
        self.x=x
        self.y=y
        self.cellx=0
        self.celly=0
        self.vx=vx
        self.vy=vy
        self.next_node=None
    
    def set_next_node(self,n):
        self.next_node=n
    
    def get_next_node(self):
        return self.next_node
    
    def set_position(self,x,y):
        self.x=x
        self.y=y
    def get_position(self):
        return self.x, self.y
    
    def set_velocity(self,vx,vy):
        self.vx=vx
        self.vy=vy
    def get_velocity(self):
        return self.vx,self.vy

class LinkedList(object):
    def __init__(self,r=None): #,r=None
        self.root=None
        self.natoms=0
    
    def add(self,new_node):
        if new_node is not None:
            new_node.next_node = self.root
            self.root = new_node
    
    def getSize(self):
        return self.natoms

    def printList(self):
        node = self.root
        while node is not None:
            listStr="LinkedList->printList() x={} y={} cellx={} celly={}"
            print(listStr.format(node.x,node.y,node.cellx,node.celly))
            node = node.next_node

class LinkedListMat(object):
    def __init__(self):
        self.pixel_root=None
        self.nx=0       #nx elements
        self.ny=0       #ny elements
        self.L=20

    def setGrid(self,nx,ny):
        linked_list=LinkedList()
        self.pixel_root=[[linked_list for j in range(nx)] for i in range(ny)]
        self.nx=nx
        self.ny=ny
        #grid of linked list
        for i in range(nx):
            for j in range(ny):
                linked_list=LinkedList()
                self.pixel_root[i][j]=linked_list
    
    def setNatoms(self,n):
        #add n atoms
        linked_list = LinkedList()
        for nn in range(n):
            node = Node()
            node.x=random.random()*(self.L)
            node.y=random.random()*(self.L)
            cellx=int(node.x*(self.nx/(1.*self.L) ))
            celly=int(node.y*(self.ny/(1.*self.L) ))
            node.cellx=cellx
            node.celly=celly
            self.pixel_root[cellx][celly].add(node)

    def printList(self):
        print("LinkedListMat->printList()")
        #node=self.root
        for n_x in range(0,self.nx):
            for n_y in range(0,self.ny):
                linked_list = self.pixel_root[n_x][n_y]
                linked_list.printList()

def main():
    n=3
    nx=10
    ny=10
    
    linkedListMatrix=LinkedListMat()
    linkedListMatrix.setGrid(nx,ny)
    linkedListMatrix.setNatoms(n)
    linkedListMatrix.printList()
    
    '''
    LinkedListTest = LinkedList()
    #node1=Node(1,1,1,1,1)
    node2=Node(2,2,2,2,2)
    #node3=Node(3,3,3,3,3)
    LinkedListTest.add(None)
    LinkedListTest.add(node2)
    #LinkedListTest.add(node3)
    LinkedListTest.printList()'''


if __name__ == "__main__":
    main()

