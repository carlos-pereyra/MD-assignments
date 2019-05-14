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
    def __init__(self,x,y,vx,vy,n=None):
        self.x=x
        self.y=y
        self.vx=vx
        self.vy=vy
        self.next_node=n
    
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
    def __init__(self,r=None):
        self.root=r
        self.size=0
    
    def add(self,x,y,vx,vy):
        new_node=Node(x,y,vx,vy,self.root)
        self.root=new_node
        self.size+=1
    
    def get_size(self):
        return self.size

    def print_list(self):
        node=self.root
        while node:
            print("x: {} y: {}".format(node.x,node.y))
            node=node.next_node


ndx=10
ndy=10
natoms=1
linked_list_matrix = [[LinkedList() for j in range(ndx)] for i in range(ndy)]


for cell_y in range(ndy):
    for cell_x in range(ndx):
        for nparticle in range(natoms):
            print("x: {} y: {} atoms: {}".format(cell_x,cell_y,nparticle))
            linked_list_matrix[cell_x][cell_y].add(8,9,1,1)
            #linked_list_matrix.remove_particles()

print(linked_list_matrix[0][9].print_list())


'''
cell=LinkedList()
cell.add(1,2,1,1)
cell.add(1,3,1,1)
cell.add(1,4,1,1)
print("test")
cell.print_list()'''

