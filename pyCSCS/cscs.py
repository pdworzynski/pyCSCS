#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from itertools import combinations,islice
from skbio.stats.ordination import pcoa
from skbio import DistanceMatrix
from scipy.sparse import dok_matrix,csc_matrix
import skbio
import pandas as pd
import biom
from multiprocessing import Process, Queue 
from collections import defaultdict
import itertools
import scipy.sparse
import pickle
import time

def filter_matrix(matrix, threshold=0.6):
    Xindices, Yindices = matrix.nonzero()
    for i in range(len(Xindices)):
        x, y = Xindices[i], Yindices[i]
        if matrix[x, y] < threshold:
            matrix[x, y] = 0.0
    matrix.eliminate_zeros()
    

def read_css(css_edges, observationids, p, cosine_threshold):
    edgesdok = dok_matrix((p, p), dtype=np.float32)
    fh = open(css_edges, "r")
    for line in fh.readlines():
        if line.find("CLUSTERID1") > -1:
            continue
        linesplit = line.split("\t")
        if float(linesplit[4]) < cosine_threshold:
            edgesdok[int(observationids[linesplit[0]]), int(observationids[linesplit[1]])] = 0.0
        else:
            edgesdok[observationids[linesplit[0]], observationids[linesplit[1]]] = float(linesplit[4])
            edgesdok[observationids[linesplit[1]], observationids[linesplit[0]]] = float(linesplit[4])
    fh.close()
    edgesdok.setdiag(1)
    return(edgesdok)

def read_bucket(features_file, normalization):
    bucket = pd.read_csv(features_file, sep = "\t",index_col=0, header=0 )
    featureids = list(bucket.index)
    sample_names = list(bucket.columns.values)

    if normalization == True:
        bucket = bucket.div(bucket.sum(axis=0), axis=1)

    return((sample_names, featureids,csc_matrix(bucket.values)))
    
def cscs_from_files(features_file, css_edges, cosine_threshold = 0.6, normalization = True, weighted = True, cpus = 1, chunk = 2):
    """ Compute CSCS from input files in GNPS buckettable and egdes format
    Args:
        features_file (table): A path to a buckettable file from GNPS
        edges (table): The CSS matrix as a edges file from GNPS
        cosine_threshold (float): Threshold under which all entries in the CSS matri will be set to 0. Set to 0.6 by default as in Sedio et al.
        normalization (boolean): Total Ion Current Sum Normalization or not
        weight (boolean): Weight all intensities or treat them as presence/absence
        cpus (int): Number of processes to run default = 1
        chunk (int): Number of samples to process in each process

    """
    sample_names, featureids, features = read_bucket(features_file, normalization)
    observationids = {str(x):index for index, x in enumerate(featureids)}
    edgesdok = read_css(css_edges, observationids, features.shape[0], cosine_threshold)

    if weighted == False:
        features = features.pa #convert to scipy sparse

    return(cscs(features, edgesdok, sample_names, cosine_threshold = cosine_threshold, normalization = normalization, weighted = weighted, cpus = cpus, chunk = chunk))

def cscs(features, edges, sample_names, cosine_threshold = 0.6, normalization = True, weighted = True, cpus = 1, chunk = 2):
    """ Compute CSCS distance
    
    Computes a CSCS distance matrix from a features matrix and a matrix of cosine similarities

    Args:
        features (table): A pandas dataframe or numpy array of feature intensities
        edges (table): The CSS matrix as a scipy.sparse.csc_matrix
        sample_names (list): list of sample names
        cosine_threshold (float): Threshold under which all entries in the CSS matri will be set to 0. Set to 0.6 by default as in Sedio et al.
        normalization (boolean): Total Ion Current Sum Normalization or not
        weight (boolean): Weight all intensities or treat them as presence/absence
        cpus (int): Number of processes to run default = 1
        chunk (int): Number of samples to process in each process

    """
    filter_matrix(features, 0.6)
    dist = parallel_make_distance_matrix(features, edges,  sample_names, cpus, chunk)
    dist = 1 - dist
    return(dist)

def compute_sum(sampleA, sampleB, edges):
    outer = sampleA.multiply(sampleB.transpose())
    finaldok = outer.multiply(edges)
    if finaldok.nnz == 0:
        return(0)
    else:
        return(sum(finaldok.data))

def single_distance(sampleA, sampleB, edges):
    """ Compute the distance between one pair of samples
    """
    start = time.time()
    cssab = compute_sum(sampleA, sampleB, edges)
    cssaa = compute_sum(sampleA, sampleA, edges)
    cssbb = compute_sum(sampleB, sampleB, edges)
    new = time.time()
    return(cssab/max(cssaa, cssbb))

def worker(input, output, edges):
    for worker_samples in iter(input.get, "STOP"):
        results=[]
        for index, sampleA, sampleB in worker_samples:
            results.append((index, single_distance(sampleA, sampleB, edges)))
        output.put([(index, result) for index, result in results])

def split_every(n, iterable):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))

def parallel_make_distance_matrix(features, edges,  sample_names, cpus, chunk):
    #Parallel stuff
    NUMBER_OF_PROCESSES = cpus
    work_chunk_size = chunk
    
    # Create queues
    task_queue = Queue()
    done_queue = Queue()
    
    #Scientific stuff
    dist =  np.zeros([features.shape[1], features.shape[1]])
    
    feature_combinations = itertools.combinations(range(0,features.shape[1]), 2)
    comb_chunks_list = list([chunk for chunk in split_every(work_chunk_size, feature_combinations)])

    for chunk_index, chunk_comb in enumerate(comb_chunks_list):
        task_queue.put([(chunk_index*work_chunk_size + comb_index, features[:,comb[0]], features[:,comb[1]]) for comb_index, comb in enumerate(chunk_comb)])
        
    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue, done_queue, edges)).start()
        
    indexed_distances = []
    for i in range(len(comb_chunks_list)):
        indexed_distances.extend(done_queue.get())
        
    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
        
    distances = list([d for i, d in sorted(indexed_distances, key=lambda id_tuple: id_tuple[0])])
    
    xs,ys = np.triu_indices(dist.shape[0],k=1)
    dist[xs,ys] = distances
    dist[ys,xs] = distances
    dist[ np.diag_indices(dist.shape[0]) ] = 1
    distdf = pd.DataFrame(dist, sample_names, sample_names)
    return(distdf)

