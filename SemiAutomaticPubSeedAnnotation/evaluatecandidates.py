
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
import scipy.stats as st

from Multi_project_scripts import query_fasta_sequences, get_pubseed_fig_information, \
                                    get_multiple_fasta_information, get_orgid_from_fig, \
                                    get_peg_candidates_from_multithread, multithreaded, \
                                    get_fig_from_genome_feature, continuous_figs, \
                                    get_featureID_frim_fig, gene_order_conservation, \
                                    get_average_normalized_bit_score, \
                                    merge_two_dicts, \
                                    update_average_normalized_bit_score, \
                                    BBHs_order_conservation, \
                                    get_gsp_confidence_score

from pubseedfig import pubseedfeature
from identifycandidates import build_candidate_profiling_vector
from Similarity_measures import jaccard_similarity


def evaluate_candidate_cluster_Bayesian(ref_df,
                                        ref_y,
                                        df_test):

    clf = ComplementNB()
    clf.fit(ref_df, ref_y)

    # Make the prediction
    y_pred = clf.predict(df_test)
    yprob = clf.predict_proba(df_test)

    df_pred = pd.DataFrame(data=yprob, index=list(df_test.index),columns=list(clf.classes_))
    return df_pred, y_pred

def evaluate_candidate_cluster_spearson(ref_df,
                                        ref_y,
                                        df_test):

    ref_fig_list = list(ref_df.index)
    candidate_fig_list = list(df_test.index)

    pvalues = []
    for fig in ref_fig_list:
        if fig:
            ref_profiling_fig = ref_df.loc[fig, :]
            for candidate_fig in candidate_fig_list:
                candidate_profiling_fig = df_test.loc[candidate_fig, :]
                stats = st.ranksums(ref_profiling_fig, candidate_profiling_fig)
                stats = st.kruskal(ref_profiling_fig, candidate_profiling_fig)
                print(fig, candidate_fig, stats.pvalue)
                pvalues.append(stats.pvalue)
    pvalues = list(set(pvalues))
    pprint(pvalues)
    print(min(pvalues))

def evaluate_candidate_cluster_decision_tree(ref_df,
                                             ref_y,
                                             df_test):
    from sklearn import tree

    clf = tree.DecisionTreeClassifier()
    clf.fit(ref_df, ref_y)

    # Make the prediction
    y_pred = clf.predict(df_test)
    yprob = clf.predict_proba(df_test)

    df_pred = pd.DataFrame(data=yprob, index=list(df_test.index),columns=list(clf.classes_))
    return df_pred, y_pred

def evaluate_candidate_cluster_(ref_df,
                                ref_y,
                                df_test):

    ref_fig_list = list(ref_df.index)
    candidate_fig_list = list(df_test.index)

    print()

def format_profiling_vector(df):

    variable_2remove = ['peg', 'operon', 'cluster', 'functional_role', 'strand']
    df = df.loc[:, [elt for elt in list(df) if elt not in variable_2remove]]
    df.ev = df.ev.apply(lambda x: True if x else False)
    df = df.set_index('ID')
    return df

def update_candidate_PEGs(candidate_fig_object_dict,
                          ref_fig_object_dict):

    # Merge both dictionaries
    fig_object_dict = merge_two_dicts(candidate_fig_object_dict,
                                      ref_fig_object_dict)
    fig_object_dict = update_average_normalized_bit_score(fig_object_dict)
    candidate_fig_object_dict = {k: v for k,v in fig_object_dict.items() if k in candidate_fig_object_dict}
    return candidate_fig_object_dict

def evaluate_candidate_clusters(ref_df,
                                candidate_cluster_fig_dict):
    '''
    :param operon_clusters_content:
    :param candidate_clusters_items:
    :return:
    '''
    ref_y = ref_df.peg

    assert type(ref_df) == pd.core.frame.DataFrame
    assert type(ref_y) == pd.core.series.Series

    ref_df = format_profiling_vector(ref_df)

    predicted_fig = {}
    for candidate_cluster in candidate_cluster_fig_dict:
        if candidate_cluster =='fig|339860.6.peg.122_124':

            pprint(candidate_cluster_fig_dict[candidate_cluster])
            '''generate_profiling_candidate_cluster'''

            candidate_df, variables = build_candidate_profiling_vector(candidate_cluster_fig_dict[candidate_cluster], ref_df)
            candidate_df = format_profiling_vector(candidate_df)

            evaluate_candidate_cluster_(ref_df, ref_y, candidate_df)
            variables = pickle.load(open('Inputs/variables.pkl', "rb"))

            '''Bayesian prediction'''
            df_pred, y_pred = evaluate_candidate_cluster_Bayesian(ref_df, ref_y, candidate_df)
            df_pred, y_pred = evaluate_candidate_cluster_decision_tree(ref_df, ref_y, candidate_df)
            print(df_pred)
            print(y_pred)

            evaluate_candidate_cluster_spearson(ref_df, ref_y, candidate_df)
            evaluate_candidate_cluster_decision_tree(ref_df, ref_y, candidate_df)
            df_pred, y_pred = evaluate_candidate_cluster_decision_tree(ref_df, ref_y, candidate_df)

            fig_list = list(df_pred.index)
            for i in range(len(fig_list)):
                fig = fig_list[i]
                predicted_role = y_pred[i]
                print(df_pred[predicted_role][fig])
                if fig not in predicted_fig:
                    predicted_fig[fig] = predicted_role

    return predicted_fig

def get_cluster_functional_roles(candidate_fig_object_dict,
                                 ref_fig_object_dict,
                                 candidate_cluster_match_dict):

    for candidate_cluster in candidate_cluster_match_dict:
        similarity_candidate_cluster = candidate_cluster_match_dict[candidate_cluster]
        matching_ref_clusters = {k:v for k,v in similarity_candidate_cluster.items() if v > 0.8}

        candidate_functional_roles = []

        for candidate_fig in candidate_fig_object_dict[candidate_cluster]:
            candiddate_fig_object = candidate_fig_object_dict[candidate_cluster][candidate_fig]
            for ref_cluster in matching_ref_clusters:
                for ref_fig in ref_fig_object_dict[ref_cluster]:
                    if ref_fig in candiddate_fig_object.peg_candidates:
                        if candiddate_fig_object.peg_candidates[ref_fig]['orthologuous']:
                            print(candidate_fig, ref_fig)
                            functional_role = ref_fig_object_dict[ref_cluster][ref_fig].functional_role
                            print(functional_role)

        # annotation_confidence_score
        print()

def get_cluster_average_normalized_bit_score(fig_object_dict1,
                                             fig_object_dict2):

    score = sum([fig_object_dict1[fig].peg_candidates[query]['average_normalized_bit_score']
                  for fig in fig_object_dict1
                  for query in fig_object_dict1[fig].peg_candidates
                  if query in fig_object_dict2])
    return score

def get_cluster_fig_dict(fig_object_dict):

    cluster_fig_dict = {}
    for fig, fig_object in fig_object_dict.items():
        if fig_object.cluster not in cluster_fig_dict:
            cluster_fig_dict[fig_object.cluster] = []
        cluster_fig_dict[fig_object.cluster].append(fig)
    return cluster_fig_dict

def get_ref_candidate_cluster_combinations(candidate_fig_object_dict,
                                           ref_fig_object_dict):

    # get candidate_cluster
    candidate_cluster_fig_dict = get_cluster_fig_dict(candidate_fig_object_dict)
    ref_cluster_fig_dict = get_cluster_fig_dict(ref_fig_object_dict)

    # Find ref clusters corresponding to candidate_clusters
    candidate_cluster_ref = {}
    for candidate_cluster in candidate_cluster_fig_dict:
        if candidate_cluster not in candidate_cluster_ref:
            candidate_cluster_ref[candidate_cluster] = []

        for fig_candidate in candidate_cluster_fig_dict[candidate_cluster]:
            for query in candidate_fig_object_dict[fig_candidate].peg_candidates:
                if query in ref_fig_object_dict and ref_fig_object_dict[query].cluster != candidate_cluster:
                    candidate_cluster_ref[candidate_cluster].append(ref_fig_object_dict[query].cluster)

    return candidate_cluster_fig_dict, ref_cluster_fig_dict, candidate_cluster_ref

def evaluate_candidate_ref_similarity(fig_object_dict1,
                                      fig_object_dict2,
                                      fig_object_dict,
                                      orthologuous_thresholds):


    # TODO(incorporate average_similarity_score)
    # TODO(fig_object_dict1 contains average_similarit_score)


    orthologues_order_conservation_score, \
    contiguous_orthologues_number,\
    orthologue_number, \
    orthologues = gene_order_conservation(fig_object_dict1,
                                          fig_object_dict2,
                                          orthologuous_thresholds)

    BBHs_order_conservation_score, \
    contiguous_BBHs_number, \
    BBHs_number, \
    BBHs = BBHs_order_conservation(fig_object_dict1,
                                   fig_object_dict2,
                                   fig_object_dict,
                                   orthologuous_thresholds)

    if orthologues_order_conservation_score or BBHs_order_conservation_score:
        length_average = (len(fig_object_dict2)
                                  + len(fig_object_dict1)) / 2
        length_diff = 1 + abs(len(fig_object_dict2)
                                  - len(fig_object_dict1))
        cluster_length_coverage = length_average / length_diff

        orthologues_order_conservation_score = orthologues_order_conservation_score * cluster_length_coverage
        BBHs_order_conservation_score = BBHs_order_conservation_score * cluster_length_coverage

        ortho_BBHs_combinaisons = list(itertools.product(orthologues, BBHs))
        ortho_BBHs = [comb[0] for comb in ortho_BBHs_combinaisons if comb[0] == comb[1]]
        ortho_BBHs_score = len(ortho_BBHs)

        gsps, gsp_confidence_score = get_gsp_confidence_score(fig_object_dict2)

        sdf = [fig_object_dict1[fig].peg_candidates[query]['average_normalized_bit_score']
               for fig in fig_object_dict1
               for query in fig_object_dict1[fig].peg_candidates
               if query in fig_object_dict2
               ]
        average_normalized_bit_score_clusters = sum(sdf)

        in_array = [orthologues_order_conservation_score,
                    BBHs_order_conservation_score,
                    gsp_confidence_score,
                    ortho_BBHs_score,
                    average_normalized_bit_score_clusters
                    ]

        out_array = map(lambda x: np.exp(x), in_array)
        out_array = [elt for elt in out_array]
        cluster_order_conservation_score = np.prod(out_array)
        return cluster_order_conservation_score, orthologues, BBHs, ortho_BBHs
    return 0.0, orthologues, BBHs, []

def get_candidate_all_ref_similarities(fig_object_dict1,
                                       fig_object_dict,
                                       ref_cluster_list,
                                       ref_fig_object_dict,
                                       ref_cluster_fig_dict,
                                       orthologuous_thresholds):

    '''
    Aim: Build a dictionary describing the similarity between cluster1 as fig_object_dict1 and
    all reference clusters in ref_fig_object_dict that match the cluster1 in ref_cluster_list
    :param fig_object_dict1:
    :param fig_object_dict:
    :param ref_cluster_list:
    :param ref_fig_object_dict:
    :param ref_cluster_fig_dict:
    :param orthologuous_thresholds:
    :return:
    '''

    ref_cluster_info = {}
    for ref_cluster in ref_cluster_list:
        fig_object_dict2 = {k: v for k, v in ref_fig_object_dict.items()
                            if k in ref_cluster_fig_dict[ref_cluster]
                            }

        cluster_order_conservation_score, \
        orthologues, \
        BBHs, \
        ortho_BBHs = evaluate_candidate_ref_similarity(fig_object_dict1,
                                                        fig_object_dict2,
                                                        fig_object_dict,
                                                        orthologuous_thresholds)

        if cluster_order_conservation_score and orthologues:
            # If there is at least 1 orthologue between two clusters

            out_dict = {'cluster_order_conservation_score': cluster_order_conservation_score,
                        'orthologue_number': len(orthologues),
                        'BBHs_number': len(BBHs),
                        'orthologues': orthologues,
                        'BBHs': BBHs,
                        'ortho_BBHs': ortho_BBHs
                        }
            ref_cluster_info[ref_cluster] = out_dict

    return ref_cluster_info

def get_proba_operon(fig_object_dict1,
                     ref_fig_object_dict,
                     ref_cluster_dict
                     ):

    m = len(fig_object_dict1)
    n = len(ref_fig_object_dict.values())
    proba = m*n*2**(-ref_cluster_dict['cluster_order_conservation_score'])
    return proba

def evaluate_candidate_fig(fig,
                           fig_object_dict1,
                           filtered_ref_cluster_fig_list,
                           ref_fig_object_dict
                           ):

    fig_matches = set(list(fig_object_dict1[fig].peg_candidates.keys())).intersection(filtered_ref_cluster_fig_list)
    orthologous = [fig2 for fig2 in fig_matches if fig_object_dict1[fig].peg_candidates[fig2]['orthologuous']]
    peg_matches = [ref_fig_object_dict[fig2].peg.title() for fig2 in orthologous]

    # Check the peg_matches consistency
    peg_matches_consistency = Counter(peg_matches)

    if len(peg_matches_consistency) == 1:
        return list(peg_matches_consistency.keys())[0]

def evaluate_candidate_cluster(fig_object_dict1,
                               fig_object_dict,
                               ref_cluster_list,
                               ref_fig_object_dict,
                               ref_cluster_fig_dict,
                               thresholds):

    '''
    Aim : The candidate_cluster is compared to all reference_clusters.
    (1) Different similarity scores are calculated to evaluate the similarity
    (2) Sort reference_clusters by score
    (3) Select highest scores and return them if score > threshold

    :param fig_object_dict1:
    :param fig_object_dict:
    :param ref_cluster_list:
    :param ref_fig_object_dict:
    :param ref_cluster_fig_dict:
    :param orthologuous_thresholds:
    :return:
    '''

    # Build a dictionary of different variables describing the comparison between rrference_cluster and the
    # candidate_cluster (BBHs, orthologues, gene order consistency)
    ref_cluster_info = get_candidate_all_ref_similarities(fig_object_dict1,
                                                          fig_object_dict,
                                                          ref_cluster_list,
                                                          ref_fig_object_dict,
                                                          ref_cluster_fig_dict,
                                                          thresholds['orthologuous_thresholds'])
    if ref_cluster_info:
        print('ref_cluster_list')
        pprint(ref_cluster_list)
        print('ref_cluster_info')
        pprint(ref_cluster_info)
        print('--------------------------------')

        # TODO(check the assignment of the candidate fig is correct)
        # TODO(determine why not better homologs)

        # TODO(check that fig|188937.1.peg.3046 and its cluster and within fig|518636.5.peg.1034 cluster)
        # 188937.1 fig|188937.1.peg.3046
        for fig in fig_object_dict1:
            for query in fig_object_dict1[fig].peg_candidates:
                if query in ref_cluster_fig_dict['fig|188937.1.peg.3045_3047']:
                    pprint(fig_object_dict1[fig].peg_candidates[query])

        # from ref_cluster_info, get the cluster_order_conservation_score
        ref_cluster_proba = {ref_cluster: get_proba_operon(fig_object_dict1,
                                                           ref_fig_object_dict,
                                                           ref_dict)
                             for ref_cluster, ref_dict in ref_cluster_info.items()
                             }

        # Reduce low Pvalue by the evalue threshold
        min_value = min(ref_cluster_proba.values())
        filtered_ref_cluster_proba = {ref_cluster: proba
                                      for ref_cluster, proba in ref_cluster_proba.items()
                                      if proba <= min_value / thresholds['prot_withdrawal']
                                      }

        filtered_ref_cluster_fig_list = [fig for fig, fig_object in ref_fig_object_dict.items()
                                         if fig_object.cluster in filtered_ref_cluster_proba]

        # for each fig in canidate_fig cluster, get all orthologous belonging to filtered_ref_cluster
        pprint(filtered_ref_cluster_proba)
        print(filtered_ref_cluster_fig_list)
        fig_assignment_dict = {}
        for fig in fig_object_dict1:
            fig_assignment = evaluate_candidate_fig(fig,
                                                    fig_object_dict1,
                                                    filtered_ref_cluster_fig_list,
                                                    ref_fig_object_dict
                                                    )
            print(fig, fig_assignment)
            fig_assignment_dict[fig] = fig_assignment
        return fig_assignment_dict

def assign_candidate_with_annotations_of_ref_cluster(fig_object_dict1,
                                                     ref_fig_object_dict,
                                                     ref_cluster_info,
                                                     ref_cluster_score):

    ref_cluster, score = ref_cluster_score

    candidate_cluster_annot = {}
    for fig, fig_object in fig_object_dict1.items():
        for comb in ref_cluster_info[ref_cluster]['orthologues']:
            candidate_fig, ortho = comb
            if fig == candidate_fig:
                if fig not in candidate_cluster_annot:
                    candidate_cluster_annot[fig] = {}

                peg = ref_fig_object_dict[ortho].peg
                fig_object_dict1[fig].peg = peg
                if ortho not in candidate_cluster_annot[fig]:
                    candidate_cluster_annot[fig][ortho] = peg

    return fig_object_dict1, candidate_cluster_annot

def evaluate_candidate_cluster_list(candidate_fig_object_dict,
                                    ref_fig_object_dict,
                                    orthologuous):

    '''
    Aim: Evaluate each candidate_cluster for its similarity with reference clusters
     to determine the function of the cluster
    Process:
    :param candidate_fig_object_dict:
    :param ref_fig_object_dict:
    :param orthologuous:
    :return:
    '''

    # Update the candidate_fig_object_dict with the average_normalized_bit_score
    # candidate_fig_object_dict = update_candidate_PEGs(candidate_fig_object_dict, ref_fig_object_dict)
    fig_object_dict = merge_two_dicts(candidate_fig_object_dict, ref_fig_object_dict)

    # get all combinations between candidate and reference clusters
    candidate_cluster_fig_dict, \
    ref_cluster_fig_dict, \
    candidate_cluster_ref = get_ref_candidate_cluster_combinations(candidate_fig_object_dict,
                                                                   ref_fig_object_dict)

    # For each candidate_cluster, matching reference_clusters are evaluated for their similarity with
    # the candidate_cluster
    all_candidate_annot = {}
    candidate_cluster_ref = {k: v
                             for k,v in candidate_cluster_ref.items()
                             if k == 'fig|518636.5.peg.1033_1036'
                             }
    for candidate_cluster in candidate_cluster_ref:
        print('candidate_cluster', candidate_cluster)
        fig_object_dict1 = {k: v for k,v in candidate_fig_object_dict.items()
                            if k in candidate_cluster_fig_dict[candidate_cluster]
                            if k in candidate_cluster_fig_dict[candidate_cluster]
                            }

        ref_cluster_list = set(candidate_cluster_ref[candidate_cluster])
        fig_assignment_dict = evaluate_candidate_cluster(fig_object_dict1,
                                                         fig_object_dict,
                                                         ref_cluster_list,
                                                         ref_fig_object_dict,
                                                         ref_cluster_fig_dict,
                                                         orthologuous)

        if fig_assignment_dict:
            all_candidate_annot.update(fig_assignment_dict)
    return candidate_fig_object_dict, all_candidate_annot

if __name__ == "__main__":
    from time import time
    tstart = time()
    print("Start")

    orthologuous_thresholds = {'single': {'aa_identity': 0.4, 'length': 0.7},
                               'cluster': {'aa_identity': 0.3, 'length': 0.6}
                               }

    if True:
        ref_fig_object_dict = pickle.load(open('Inputs/pickle/ref_fig_network1_object_dict.pkl', "rb"))
        candidate_fig_object_dict = pickle.load(open('Inputs/pickle/candidate_fig_object_dict.pkl', "rb"))

        pprint(candidate_fig_object_dict)

        identify_candidates_in_genomes(genome_ID,
                                       ref_fig_object_dict,
                                       prot_list,
                                       Blast_parameters,
                                       workers,
                                       thresholds)

    else:
        print()

    tend = time()
    print('Process duration : ', int(tend - tstart))
    print ("Done")