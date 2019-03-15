
import csv
import itertools
import operator
import os
import pickle
import urllib
from collections import Counter
from html.parser import HTMLParser
from itertools import groupby
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from os import startfile
from pprint import pprint
from random import seed, randint, SystemRandom, sample
from re import findall, match
from subprocess import Popen
from time import sleep
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from webbrowser import open as webbrowser_open
from sklearn.naive_bayes import GaussianNB
from sklearn.naive_bayes import ComplementNB
# from sklearn.model_selection import train_test_split

import numpy as np
import pandas as pd

from Multi_project_scripts import query_fasta_sequences, get_pubseed_fig_information, \
                                    get_multiple_fasta_information, get_orgid_from_fig, \
                                    get_peg_candidates_from_multithread, multithreaded, \
                                    get_fig_from_genome_feature, continuous_figs, \
                                    get_featureID_frim_fig, get_fig_clusters, get_candidate_clusters_for_operon, \
                                    get_multithreading_protein_significant_similarities, \
                                    merge_two_dicts, \
                                    update_average_normalized_bit_score, \
                                    retreive_fig_self_similarity, \
                                    get_sequence_similarity_genome_with_refPEGs, \
                                    update_refPEGs, \
                                    get_genomes_fig_to_match, \
                                    get_pubseed_single_fig_information

from pubseedfig import pubseedfeature
from decorators import function_timing

@function_timing
def get_potential_candidate_clusters(genome_id,
                                     fig_object_dict,
                                     my_pegs_interest,
                                     values,
                                     thresholds,
                                     workers):
    '''
    Aim : Determine if the genome contains potential candidate proteins of my_pegs_interest
    :param genome_id:
    :param fig_object_dict:
    :param my_pegs_interest:
    :param values:
    :param thresholds:
    :param workers:
    :return: True / False
    '''

    fig_in_peg_interest = [fig for fig in fig_object_dict if fig_object_dict[fig].peg in my_pegs_interest]
    fig_sequences = {fig_object_dict[fig].ID: [fig_object_dict[fig].prot_sequence] for fig in fig_in_peg_interest}

    query_fig_candidates = multithreaded(genome_id, fig_sequences, values, workers)
    evalue_list = [value2[0] for value in query_fig_candidates.values() for value2 in value.values()]
    return any(x <= thresholds['prot_withdrawal'] for x in evalue_list)

@function_timing
def get_candidate_cluster(genome_id, fig_object_dict, Blast_parameters, thresholds, workers):

    reference_cluster_list = set(list([fig_object_dict[fig].cluster for fig in fig_object_dict]))

    for reference_cluster in reference_cluster_list:
        if reference_cluster:
            fig_list = set(list([fig for fig in fig_object_dict if fig_object_dict[fig].cluster == reference_cluster]))
            fig_sequences = {fig_object_dict[fig].ID: [fig_object_dict[fig].prot_sequence] for fig in fig_list}

            query_fig_candidates = get_multithreading_protein_significant_similarities(genome_id,
                                                                                       fig_sequences,
                                                                                       Blast_parameters,
                                                                                       workers)

            candidate_fig_list = list(set([k1 for k, v in query_fig_candidates.items() for k1 in list(v.keys()) if v]))
            fig_query_evalues = {fig: {query: v[fig] for query, v in query_fig_candidates.items() if fig in v.keys()} for fig in
                                 candidate_fig_list}

            # Find candidates clusters from fig_candidates
            candidate_clusters = get_candidate_clusters_for_operon(query_fig_candidates, thresholds)

            # built similarity between reference operon and candidate clusters : add to candidate cluster the fig_candidates result given fig content for the cluster
            for candidate_cluster in candidate_clusters:
                end_feature_ID = get_featureID_frim_fig(candidate_cluster[-1])
                candidate_cluster_id = '{0}_{1}'.format(candidate_cluster[0], end_feature_ID)

                if candidate_cluster_id not in candidate_clusters_items:
                    candidate_clusters_items[candidate_cluster_id] = {}

                for fig in candidate_cluster:
                    if fig in fig_query_evalues:
                        if fig not in candidate_clusters_items[candidate_cluster_id]:
                            candidate_clusters_items[candidate_cluster_id][fig] = {}
                        query = list(fig_query_evalues[fig].keys())[0]
                        candidate_clusters_items[candidate_cluster_id][fig][query] = fig_query_evalues[fig][query]

    return candidate_clusters_items

@function_timing
def identify_potential_clusters(genome_id,
                                fig_object_dict,
                                Blast_parameters,
                                thresholds,
                                workers):

    ref_fig_list = list(fig_object_dict.keys())
    ref_fig_sequences = {fig_object_dict[fig].ID: [fig_object_dict[fig].prot_sequence] for fig in ref_fig_list}

    query_fig_candidates = get_multithreading_protein_significant_similarities(genome_id,
                                                                               ref_fig_sequences,
                                                                               Blast_parameters,
                                                                               workers)

    # Find candidates clusters from fig_candidates
    candidate_clusters = get_candidate_clusters_for_operon(query_fig_candidates,
                                                           thresholds)
    candidate_fig_list = list(set([k1 for k, v in query_fig_candidates.items()
                                   for k1 in list(v.keys())
                                   if v]))
    fig_query_evalues = {fig: {query: v[fig] for query, v in query_fig_candidates.items()
                               if fig in v.keys()}
                         for fig in candidate_fig_list}

    candidate_clusters_items = {}
    # built similarity between reference operon and candidate clusters : add to candidate cluster the fig_candidates result given fig content for the cluster
    for candidate_cluster in candidate_clusters:
        end_feature_ID = get_featureID_frim_fig(candidate_cluster[-1])
        candidate_cluster_id = '{0}_{1}'.format(candidate_cluster[0], end_feature_ID)

        if candidate_cluster_id not in candidate_clusters_items:
            candidate_clusters_items[candidate_cluster_id] = {}

        for fig in candidate_cluster:
            if fig in fig_query_evalues:
                if fig not in candidate_clusters_items[candidate_cluster_id]:
                    candidate_clusters_items[candidate_cluster_id][fig] = {}
                for query in list(fig_query_evalues[fig].keys()):
                # query = list(fig_query_evalues[fig].keys())[0]
                    candidate_clusters_items[candidate_cluster_id][fig][query] = fig_query_evalues[fig][query]

    return candidate_clusters_items

@function_timing
def filter_candidate_clusters(candidate_clusters_items,
                              thresholds):

    '''
    Aim: For each candidate_cluster, reference_clusters found are discarded
    if there is not a single reference_fig in cluster whose similarity is above thresholds['filter']

    :param candidate_clusters_items:
    :param thresholds:
    :return:
    '''

    fig_list = list(set([fig for v in candidate_clusters_items.values() for fig in v.keys()]))
    cluster_groups = get_fig_clusters(fig_list, 5)
    cluster_groups = [group for group in cluster_groups]

    for candidate_cluster in list(candidate_clusters_items.keys()):
        if candidate_cluster:
            candidate_cluster_fig_list = [fig for fig in candidate_clusters_items[candidate_cluster]]
            # Remove duplicate clusters with lesser fig content

            for group in cluster_groups:
                inter = set(candidate_cluster_fig_list).intersection(group)
                if inter and len(inter)<len(group):
                    del candidate_clusters_items[candidate_cluster]

        # In the cluster, if any similarity is relevant, the cluster is kept for further analysis
        if candidate_cluster in list(candidate_clusters_items.keys()):
            cluster_relevancy = False
            for fig in candidate_cluster_fig_list:
                if not cluster_relevancy:
                    e_values = [scores['evalue']
                                for query,scores in candidate_clusters_items[candidate_cluster][fig].items()]
                    if min(e_values) <= thresholds['filter']:
                        cluster_relevancy = True
            if not cluster_relevancy:
                del candidate_clusters_items[candidate_cluster]
    return candidate_clusters_items

def instantiate_candidate_fig_object(fig, candidate_cluster, peg_candidates):

    operon = ''
    peg = ''
    fig_seq_info = query_fasta_sequences(fig)
    fig_information = get_pubseed_single_fig_information(fig)

    fig_object = pubseedfeature(fig,
                                fig_information['Function'],
                                fig_information['Size'],
                                fig_information['Start'],
                                fig_information['End'],
                                peg_candidates,
                                operon,
                                candidate_cluster,
                                peg,
                                fig_seq_info['Seq'],
                                fig_information['Gap'],
                                fig_information['Ev'],
                                fig_information['strand'])
    return fig, candidate_cluster, fig_object

@function_timing
def instantiate_multiple_candidate_fig(candidate_clusters_items, workers):

    fig_list = []
    candidate_cluster_list = []
    peg_candidates_list = []
    for candidate_cluster in candidate_clusters_items:
        for fig in candidate_clusters_items[candidate_cluster]:
            fig_list.append(fig)
            candidate_cluster_list.append(candidate_cluster)
            peg_candidates_list.append(candidate_clusters_items[candidate_cluster][fig])

    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results

    args = zip(fig_list,
               candidate_cluster_list,
               peg_candidates_list)

    results = pool.starmap(instantiate_candidate_fig_object, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    candidate_cluster_fig_dict = {}
    for fig_elements in results:
        fig, candidate_cluster, fig_object = fig_elements
        if candidate_cluster not in candidate_cluster_fig_dict:
            candidate_cluster_fig_dict[candidate_cluster] = {}
        candidate_cluster_fig_dict[candidate_cluster][fig] = fig_object
    return candidate_cluster_fig_dict

@function_timing
def instantiate_candidate_fig(candidate_clusters_items):

    candidate_cluster_fig_dict = {}
    for candidate_cluster in candidate_clusters_items:
        for fig in candidate_clusters_items[candidate_cluster]:
            operon = ''
            peg = ''
            fig_seq_info = query_fasta_sequences(fig)
            all_fig_information = get_pubseed_single_fig_information(fig)
            fig_information = all_fig_information[fig]
            peg_candidates = candidate_clusters_items[candidate_cluster][fig]

            fig_object = pubseedfeature(fig,
                                   fig_information['Function'],
                                   fig_information['Size'],
                                   fig_information['Start'],
                                   fig_information['End'],
                                   peg_candidates,
                                   operon,
                                   candidate_cluster,
                                   peg,
                                   fig_seq_info['Seq'],
                                   fig_information['Gap'],
                                   fig_information['Ev'],
                                   fig_information['strand'])

            if candidate_cluster not in candidate_cluster_fig_dict:
                candidate_cluster_fig_dict[candidate_cluster] = {}
            candidate_cluster_fig_dict[candidate_cluster][fig] = fig_object
    return candidate_cluster_fig_dict

def build_candidate_profiling_vector(candidate_cluster_fig, ref_df):

    variables = list(ref_df) + ['ID']
    candidate_profiling_vector = {}
    for var in variables:
        if var not in candidate_profiling_vector:
            candidate_profiling_vector[var] = []

    for fig in candidate_cluster_fig:
        fig_object = candidate_cluster_fig[fig]
        fig_variables = fig_object.__dict__
        for k,v in fig_variables.items():
            if k in variables:
                candidate_profiling_vector[k].append(v)

        similarity_variables = [elt for elt in variables if elt not in fig_variables]
        for var in similarity_variables:
            if var in fig_object.peg_candidates:
                print(fig_object.peg_candidates[var])
        print()
            # value = fig_object.peg_candidates[var][1] if var in fig_object.peg_candidates.keys() else 0.0
            # candidate_profiling_vector[var].append(value)

    candidate_df = pd.DataFrame(candidate_profiling_vector)
    return candidate_df, variables

@function_timing
def identify_candidate_clusters(genome_ID,
                                fig_object_dict,
                                my_pegs_interest,
                                Blast_parameters,
                                workers,
                                thresholds):
    '''
    :param genome_id:
    :param operon_clusters_content:
    :return:

    1. pre-screening : only screen proteins in my_pegs_interest to speed-up process
    2. If potential candidates, screen all proteins to build similarity network
        for each cluster:
            1: Blast all fig of the cluster:
            2: find clusters
            3: for each candidate cluster:
                1: instanciate a object
    '''
    # Blast only ref_fig from my_pegs_interest to estimate if genome has potential candidate of these my_pegs_interest
    has_potential_candidate = get_potential_candidate_clusters(genome_ID,
                                                               fig_object_dict,
                                                               my_pegs_interest,
                                                               Blast_parameters,
                                                               thresholds,
                                                               workers)

    if has_potential_candidate:
        # The geonme has potential candidates. Screen all proteins to build similarity network
        # Identify all potential clusters correpsonding either to my_pegs_interest or their homologs
        candidate_clusters_items = identify_potential_clusters(genome_ID,
                                                               fig_object_dict,
                                                               Blast_parameters,
                                                               thresholds,
                                                               workers)

        # Discard matches that are not relevant (low P-value)
        candidate_clusters_items = filter_candidate_clusters(candidate_clusters_items,
                                                             thresholds)

        # For all candidate_fig, instantiate a fig_object
        candidate_cluster_fig_dict = instantiate_multiple_candidate_fig(candidate_clusters_items,
                                                                        workers)
        return candidate_cluster_fig_dict

def update_cluster_information(updating_candidate_cluster):

    fig_list = list(updating_candidate_cluster.keys())
    all_fig_information = build_fig_information(fig_list[0])
    fig_list = fig_list
    for fig in fig_list:
        if fig in all_fig_information:
            fig_o = updating_candidate_cluster[fig]
            fig_o.gap = all_fig_information[fig]['Gap'] if all_fig_information[fig]['Gap'].isdigit() else ''
            fig_o.fig_start = all_fig_information[fig]['Start']
            fig_o.fig_end = all_fig_information[fig]['End']
    return updating_candidate_cluster

def update_cluster_dict(candidate_cluster_dict):

    for candidate_cluster in list(candidate_cluster_dict.keys()):
        updated_candidate_cluster = update_cluster_information(candidate_cluster_dict[candidate_cluster])
        candidate_cluster_dict[candidate_cluster] = updated_candidate_cluster
    return candidate_cluster_dict

def update_refPEGs_from_candidate_fig(genome_ID,
                                      figs_to_match,
                                      figs_sequences,
                                      Blast_parameters,
                                      workers):

    reffigs_candidates = get_sequence_similarity_genome_with_refPEGs(genome_ID,
                                                                     figs_to_match,
                                                                     figs_sequences,
                                                                     Blast_parameters,
                                                                     workers)
    return reffigs_candidates

@function_timing
def retreive_reffig_and_candidate_fig_similarities(ref_fig_object_dict,
                                                   candidate_fig_object_dict,
                                                   Blast_parameters,
                                                   workers):
    '''
    Aim: get similarity of candidate_fig and ref_fig
    Process: 1. Blast candidate_fig protein sequence into ref_genome to get similarity information
             2. Multithreading
    :param ref_fig_object_dict:
    :param candidate_fig_object_dict:
    :param Blast_parameters:
    :param workers:
    :return:
    '''

    genomes_fig_to_match = get_genomes_fig_to_match(ref_fig_object_dict)
    # get all candidate_fig protein sequences
    figs_sequences = {fig: [item.prot_sequence] for fig, item in candidate_fig_object_dict.items()}

    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results
    genomes_to_query = []
    fig_list_per_genome_to_match = []

    for key, value in genomes_fig_to_match.items():
        genomes_to_query.append(key)
        fig_list_per_genome_to_match.append(value)

    args = zip(genomes_to_query,
               fig_list_per_genome_to_match,
               itertools.repeat(figs_sequences),
               itertools.repeat(Blast_parameters),
               itertools.repeat(workers))

    results = pool.starmap(update_refPEGs_from_candidate_fig, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    for reffigs_candidates in results:
        ref_fig_object_dict = update_refPEGs(ref_fig_object_dict, reffigs_candidates)
    return ref_fig_object_dict

@function_timing
def retreive_reffig_with_candidate_fig(ref_fig_object_dict, candidate_fig_object_dict, Blast_parameters, workers):

    '''
    Aim : Blast candidate_fig on ref genomes to get similarity data
    :param ref_fig_object_dict:
    :param candidate_fig_object_dict:
    :return:
    '''

    genomes_fig_to_match = get_genomes_fig_to_match(ref_fig_object_dict)
    # get all candidate_fig protein sequences
    figs_sequences = {fig: [item.prot_sequence] for fig, item in candidate_fig_object_dict.items()}

    for genome_ID in genomes_fig_to_match:
        if genome_ID:
            figs_to_match = genomes_fig_to_match[genome_ID]
            reffigs_candidates = get_sequence_similarity_genome_with_refPEGs(genome_ID, figs_to_match, figs_sequences,
                                                                             Blast_parameters, workers)
            ref_fig_object_dict = update_refPEGs(ref_fig_object_dict, reffigs_candidates)
    return ref_fig_object_dict

def retreive_average_normalized_bit_score(ref_fig_object_dict, candidate_fig_object_dict):

    fig_object_dict = merge_two_dicts(candidate_fig_object_dict,
                                      ref_fig_object_dict)

    fig_object_dict = update_average_normalized_bit_score(fig_object_dict)
    updated_ref_fig_object_dict = {k: v for k,v in fig_object_dict.items()
                                   if k in ref_fig_object_dict
                                   }
    updated_candidate_fig_object_dict = {k: v for k,v in fig_object_dict.items()
                                         if k in candidate_fig_object_dict
                                         }
    return ref_fig_object_dict, candidate_fig_object_dict

@function_timing
def identify_candidates_in_genomes(genome_ID,
                                   ref_fig_object_dict,
                                   prot_list,
                                   Blast_parameters,
                                   workers,
                                   thresholds):
    '''
    Aim : (i) Identify potential candidates in the genome matching the reference_fig_list
          (ii) Retreive various similarity information
    Process: (i) blast the ref_fig corresponding to prot_list and keep if similarity is relevant
             (ii) blast all ref_fig (prot_list and their homologs) to get network of similarities and other info
             (iii) blast back candidate_fig to ref_genomes to retreive their similarity (required for BBHs)
    :param genome_ID:
    :param ref_fig_object_dict:
    :param prot_list:
    :param Blast_parameters:
    :param workers:
    :param thresholds:
    :return:
    '''

    # Identify candidate clusters by (i) blasting reference_fig into the genome (ii) retreive information
    # candidate_cluster_fig_dict = identify_candidate_clusters(genome_ID,
    #                                                          ref_fig_object_dict,
    #                                                          prot_list,
    #                                                          Blast_parameters,
    #                                                          workers,
    #                                                          thresholds)
    #
    # pickle.dump(candidate_cluster_fig_dict,
    #             open('Inputs/pickle/candidate_cluster_fig_dict.pkl', "wb"),
    #             protocol=pickle.HIGHEST_PROTOCOL)
    candidate_cluster_fig_dict = pickle.load(open('Inputs/pickle/candidate_cluster_fig_dict.pkl', "rb"))
    print('candidate_cluster_fig_dict is load in identify_candidates_in_genomes function')
    print('-------------------------')

    candidate_fig_object_dict = {}
    if candidate_cluster_fig_dict:
        # Reduce the level of the candidate_cluster_fig_dict :
        # From get candidate_fig_object_dict[cluster][fig] = fig_object, get candidate_fig_object_dict[fig] = fig_object
        candidate_fig_object_dict = {fig: fig_object
                                     for cluster in candidate_cluster_fig_dict
                                     for fig, fig_object in candidate_cluster_fig_dict[cluster].items()
                                     }

        # Get the similarity score of the PEG itself
        # The similarity score of the PEG itself is required for further calculation
        candidate_fig_object_dict = retreive_fig_self_similarity(genome_ID,
                                                                 candidate_fig_object_dict,
                                                                 Blast_parameters,
                                                                 workers)

        # Blast candidate_fig protein sequence into ref_genome to get similarity information
        ref_fig_object_dict = retreive_reffig_and_candidate_fig_similarities(ref_fig_object_dict,
                                                                             candidate_fig_object_dict,
                                                                             Blast_parameters,
                                                                             workers)

        # for each similarity information, calculate average_normalized_bit_score
        ref_fig_object_dict, \
        candidate_fig_object_dict = retreive_average_normalized_bit_score(ref_fig_object_dict,
                                                                          candidate_fig_object_dict)
    return ref_fig_object_dict, candidate_fig_object_dict

if __name__ == "__main__":
    from time import time
    tstart = time()
    print("Start")

    cutoff = 1e-05
    sequence = ""
    Blast_parameters = {"seq_type":"aa",
              "filter" : "F",
              "evalue" : str(cutoff),
              "wsize" : "0",
              "fasta" : sequence,
              "organism" : '1005435.3',
              "act" : "BLAST"
            }
    workers = 50

    thresholds = {'cutoff':1e-05,
                  'filter':1e-25,
                  'operon_reannotation':1e-35,
                  'prot_reannotation':1e-275,
                  'prot_withdrawal':1e-30,
                  'cluster_threshold':5
        }
    Blast_parameters['cutoff'] = thresholds['cutoff']
    genome_ID = '339860.6'

    prot_list = ['HdrC', 'MvhD', 'HdrB', 'HdrA', 'MetF', 'metV']
    prot_list = map(elt.title(), prot_list)
    print(prot_list)
    if prot_list:
        ref_fig_object_dict = pickle.load(open('Inputs/pickle/ref_fig_network_object_dict.pkl', "rb"))

        ref_fig_object_dict, candidate_fig_object_dict = identify_candidates_in_genomes(genome_ID,
                                                                                        ref_fig_object_dict,
                                                                                        prot_list,
                                                                                        Blast_parameters,
                                                                                        workers,
                                                                                        thresholds)

        # pickle.dump(ref_fig_object_dict, open('Inputs/pickle/ref_fig_network1_object_dict.pkl', "wb"),
        #             protocol=pickle.HIGHEST_PROTOCOL)
        # pickle.dump(candidate_fig_object_dict, open('Inputs/pickle/candidate_fig_object_dict.pkl', "wb"),
        #             protocol=pickle.HIGHEST_PROTOCOL)

    tend = time()
    print('Process duration : ',int(tend - tstart))
    print ("Done")
