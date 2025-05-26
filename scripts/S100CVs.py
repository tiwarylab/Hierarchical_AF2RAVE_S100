#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: xgu
"""
import numpy as np
import mdtraj as md
from collections import deque
import random
import os
#resids_anchor=np.loadtxt('ends.index')

oxy_EF1=['19-O','27-O','32-OE1','32-OE2']  ## resSeq-name
oxy_EF2=['62-OD1','66-OD2','68-O','73-OE1','73-OE2'] ## resSeq-name

def ion_oxy_pairs(top):   # for PDB files generate from rAF2
    pairs = []
    ca1 = top.select('resSeq 1 and element Ca and chainid 2')[0]
    for oxy in oxy_EF1:
        lls = oxy.split('-')
        oid = top.select(f'resSeq {lls[0]} and name {lls[1]} and chainid 0')[0]
        pairs.append([ca1, oid])

    ca1 = top.select('resSeq 2 and element Ca and chainid 2')[0]
    for oxy in oxy_EF2:
        lls = oxy.split('-')
        oid = top.select(f'resSeq {lls[0]} and name {lls[1]} and chainid 0')[0]
        pairs.append([ca1, oid])

    ca1 = top.select('resSeq 1 and element Ca and chainid 3')[0]
    for oxy in oxy_EF1:
        lls = oxy.split('-')
        oid = top.select(f'resSeq {lls[0]} and name {lls[1]} and chainid 1')[0]
        pairs.append([ca1, oid])

    ca1 = top.select('resSeq 2 and element Ca and chainid 3')[0]
    for oxy in oxy_EF2:
        lls = oxy.split('-')
        oid = top.select(f'resSeq {lls[0]} and name {lls[1]} and chainid 1')[0]
        pairs.append([ca1, oid])
    return pairs


ends = [int(i) for i in np.loadtxt('/home/xg23/scratch/S100/rave/scripts/ends.index')]
def anchor_pairDist(traj):
    anchors = traj.top.select(f'name CA and chainid 0')[ends]
    atom_pairs = []
    for i in range(len(anchors)):
        for j in range(i+1, len(anchors)):
            atom_pairs.append([anchors[i], anchors[j]])
    atom_pairs = np.asarray(atom_pairs)
    print(atom_pairs)
    distA = md.compute_distances(traj, atom_pairs, periodic=False)


    anchors = traj.top.select(f'name CA and chainid 1')[ends]
    atom_pairs = []
    for i in range(len(anchors)):
        for j in range(i+1, len(anchors)):
            atom_pairs.append([anchors[i], anchors[j]])
    atom_pairs = np.asarray(atom_pairs)
    print(atom_pairs)
    distB = md.compute_distances(traj, atom_pairs, periodic=False)
    
    return distA, distB



def closest_heavy_atoms_dist(resid1, resid2, top, traj=0):
    select0=top.select(f"resid {resid1} and sidechain and not element H")  #Residue index (0-based)
    select1=top.select(f"resid {resid2} and sidechain and not element H")  #Residue index (0-based)
    distances=md.compute_distances(traj, [[a1, a2] for a1 in select0 for a2 in select1],periodic=False)
    distances = np.min(distances, axis=1)
    return distances

def find_distal_atom(topology, residue_index):
    sidechain = topology.select(f"resid {residue_index} and sidechain and not element H") #Residue index (0-based)
    ca_atom = topology.select(f"resid {residue_index} and name CA")
        
    if len(ca_atom) == 0:
        raise ValueError(f"No alpha carbon found in residue {residue_index}")
        
    if len(sidechain) == 0:
        return ca_atom[0]
    else:        
        # Perform a BFS to find the most distal atom
        max_distance = 0
        distal_atom = ca_atom
        visited = set()
        queue = deque([(ca_atom, 0)])

        while queue:
            current_atom, distance = queue.popleft()
            if distance > max_distance:
                max_distance = distance
                distal_atom = current_atom

            for bond in topology.bonds:   
                if bond[0].index == current_atom and (bond[1].index in sidechain):
                    neighbor = bond[1].index                    
                elif bond[1].index == current_atom and (bond[0].index in sidechain):
                    neighbor = bond[0].index
                else:
                    continue

                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, distance + 1))

        return distal_atom
    
def distal_heavy_atoms_dist(resid1, resid2, top, traj=0):
    select0 = find_distal_atom(top, resid1)
    select1 = find_distal_atom(top, resid2)
    distances = md.compute_distances(traj, [[select0, select1]],periodic=False)
    return distances

id_p = np.loadtxt('/home/xg23/scratch/S100/rave/scripts/index_pocket.txt',dtype='int')
def pocket_CVs(traj):
    Lchain = len(traj.top.select('chainid 0 and name CA'))
    distA = []
    for cv in id_p:
        distA.append(np.squeeze(distal_heavy_atoms_dist(cv[0], cv[1], traj.top, traj)))
    distA = np.transpose(distA)

    distB = []
    cvs = id_p + Lchain
    for cv in cvs:
        distB.append(np.squeeze(distal_heavy_atoms_dist(cv[0], cv[1], traj.top, traj)))
    distB = np.transpose(distB)
    return distA, distB



def RegSpaceClustering(z, min_dist, max_centers=200, batch_size=100,randomseed=0,periodicity=0):
    '''Regular space clustering.
    Args:
        data: ndarray containing (n,d)-shaped float data
        max_centers: the maximum number of cluster centers to be determined, integer greater than 0 required
        min_dist: the minimal distances between cluster centers
    '''
    random.seed(5)
    num_observations, d = z.shape
    p = np.hstack((0,np.random.RandomState(seed=randomseed).permutation(num_observations-1)+1))
    data = z[p]
    center_list = data[0, :].copy().reshape(d,1)
    centerids=[p[0]+1]
    i = 1
    while i < num_observations:
        x_active = data[i:i+batch_size, :]
        differences=np.abs(np.expand_dims(center_list.T,0) - np.expand_dims(x_active,1))
        differences=np.max(np.stack((differences,periodicity-differences)),axis=0)
        distances = np.sqrt((np.square(differences)).sum(axis=-1))
        indice = tuple(np.nonzero(np.all(distances > min_dist, axis=-1))[0])
        if len(indice) > 0:
            # the first element will be used
            #print(center_list.shape,x_active.shape,x_active[indice[0]].shape)
            center_list = np.hstack((center_list, x_active[indice[0]].reshape(d,1)))
            centerids.append(p[i+indice[0]]+1)
            i += indice[0]
        else:
            i += batch_size
        if len(centerids) >= max_centers:
            print("%i centers: Exceeded the maximum number of cluster centers!\n"%len(centerids))
            print("Please increase dmin!\n")
            raise ValueError
    print("Found %i centers!"%len(centerids))
    return center_list,centerids
