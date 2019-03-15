
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
from time import time

import numpy as np
import pandas as pd
import scipy.stats as st
import gc

from Multi_project_scripts import query_fasta_sequences, get_pubseed_fig_information, \
                                    get_multiple_fasta_information, get_orgid_from_fig, \
                                    get_peg_candidates_from_multithread, multithreaded, \
                                    get_fig_from_genome_feature, continuous_figs, \
                                    get_multithreading_protein_significant_similarities, \
                                    get_average_normalized_bit_score, \
                                    get_fig_clusters, \
                                    get_featureID_frim_fig, \
                                    update_average_normalized_bit_score, \
                                    get_multiple_pubseed_fig_information, \
                                    retreive_fig_self_similarity, \
                                    get_sequence_similarity_genome_with_refPEGs, \
                                    update_refPEGs, \
                                    get_genomes_fig_to_match

from reduceref import reduce_refdata
from pubseedfig import pubseedfeature

def instantiate_reference_fig(df):

    fig_object_dict = {}
    for row in list(df.index):
        fig_information = df.loc[row,:]
        fig = df['Fig'][row]
        fig_seq_info = query_fasta_sequences(fig)
        peg_candidates = {}
        operon = ''
        fig_object = pubseedfeature(fig,
                           fig_information['Function'],
                           float(fig_information['Size'])/3,
                           fig_information['Start'],
                           fig_information['End'],
                           peg_candidates,
                           operon,
                           fig_information['Cluster'],
                           fig_information['Peg'],
                           fig_seq_info['Seq'],
                           fig_information['Gap'],
                           fig_information['Ev'],
                           fig_information['strand'])
        if fig not in fig_object_dict:
            fig_object_dict[fig] = fig_object
    return fig_object_dict

def build_seq_similarity_genome(genome_ID,
                                fig_object_dict,
                                genome_peg_fig_to_match,
                                Blast_parameters,
                                workers):

        '''
        Aim : Blast GSP sequences against the genome_ID to get similarity network
        :param genome_ID:
        :param fig_object_dict:
        :param genome_peg_fig_to_match:
        :param Blast_parameters:
        :param workers:
        :return:
        '''

        figs_to_match = genome_peg_fig_to_match
        figs_sequences = {fig: [item.prot_sequence] for fig, item in fig_object_dict.items()
                          # figID not in figs_to_match
                          }
        reffigs_candidates = get_sequence_similarity_genome_with_refPEGs(genome_ID,
                                                                         figs_to_match,
                                                                         figs_sequences,
                                                                         Blast_parameters,
                                                                         workers)
        return reffigs_candidates

def build_seq_similarity_network(fig_object_dict,
                                 Blast_parameters,
                                 workers):

    genomes_fig_to_match = get_genomes_fig_to_match(fig_object_dict)
    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results

    args = zip(list(genomes_fig_to_match.keys()),
               itertools.repeat(fig_object_dict),
               list(genomes_fig_to_match.values()),
               itertools.repeat(Blast_parameters),
               itertools.repeat(workers)
               )

    results = pool.starmap(build_seq_similarity_genome, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    for reffigs_candidates in results:
        fig_object_dict = update_refPEGs(fig_object_dict, reffigs_candidates)
    return query_candidates

def save_operon_clusters_content_2tsv(df,
                                      ref_file):

    assert type(df) == pd.core.frame.DataFrame
    df.to_csv(ref_file, sep='\t', encoding='utf-8')

def get_operon_clusters_content_from_tsv(file):

    df = pd.read_csv(file, sep='\t', encoding='utf-8')
    print(df)

def collect_pubseed_feature_information(operon_clusters_content_file,
                                        my_operon_interest):

    operon_clusters_content = set_operon_clusters_content(operon_clusters_content_file, my_operon_interest)
    return operon_clusters_content

def get_pubseed_spreadsheet(file,
                            prot_list):

    df = pd.read_csv(file, sep='\t', encoding='utf-8')

    # Reduce df on prot_list
    variables_2keep = ['Organism']
    variables_2keep.extend(prot_list)

    # Update prot_list to match variables of df
    variables_2keep = list(set(list(df)).intersection(variables_2keep))
    df = df.loc[:, variables_2keep]

    def f(x):
        if x == ' ':
            return ''
        else:
            return x

    df = df.applymap(f)
    df = df.loc[(df[prot_list] != '').any(axis=1), :]

    def g(x):
        new_x = x.split('_')[0]
        if new_x.isdigit():
            return new_x
        else:
            return x

    df = df.applymap(g)
    df['Genome_ID'] = df['Organism'].map(lambda x: findall(' \((\d+\.\d+)\)', x)[0])
    return df

def gold_standard_proteins(file):

    df = pd.read_csv(file, sep='\t', encoding='utf-8')
    return df

def get_peg_fig_df(pubseed_spreadsheet_df,
                   prot_list):

    fig_clusters = {'Peg': [], 'Fig': []}

    for row in list(pubseed_spreadsheet_df.index):
        genome_data = pubseed_spreadsheet_df.loc[row]
        feature_peg_dict = {feature:prot for prot in prot_list for feature in genome_data[prot].split(', ') if feature}
        feature_list = [feature for prot in prot_list for feature in genome_data[prot].split(', ') if feature]

        for feature in feature_list:
            fig = get_fig_from_genome_feature(genome_data['Genome_ID'], str(feature))
            if fig:
                fig_clusters['Fig'].append(fig)
                fig_clusters['Peg'].append(feature_peg_dict[feature])

    df = pd.DataFrame(fig_clusters)
    return df

def load_pubseed_spreadsheet_cluster_content(PubSeed_speadsheet_file,
                                             prot_list):
    pubseed_spreadsheet_df = get_pubseed_spreadsheet(PubSeed_speadsheet_file, prot_list)
    peg_fig_df = get_peg_fig_df(pubseed_spreadsheet_df, prot_list)
    return peg_fig_df

def merge_pubseed_spreadsheet_gsp(pubseed_spreadsheet_df,
                                  gsp_df):

    frames = [pubseed_spreadsheet_df, gsp_df]
    concat_df = pd.concat(frames)

    # Add the genome_ID information to the df
    concat_df['Genome_ID'] = concat_df.Fig.map(lambda x: get_orgid_from_fig(x))
    return concat_df

def retreive_pubseed_information(pubseed_spreadsheet_df,
                                 gsp_df,
                                 cluster_threshold,
                                 workers):

    concat_df = merge_pubseed_spreadsheet_gsp(pubseed_spreadsheet_df, gsp_df)
    concat_fig_list = list(concat_df.Fig)

    fig_dict = {'Cluster': [None]*len(concat_fig_list), 'Fig': concat_fig_list}

    for genome_ID in list(set(concat_df.Genome_ID)):
        genome_df = concat_df.loc[concat_df.Genome_ID == genome_ID]

        # For each genome, determine clusters of fig
        fig_list = list(genome_df.Fig)
        cluster_groups = get_fig_clusters(fig_list, cluster_threshold)
        for cluster in cluster_groups:
            end_feature_ID = get_featureID_frim_fig(cluster[-1])
            cluster_ID = '{0}_{1}'.format(cluster[0], end_feature_ID)

            for fig in cluster:
                fig_dict['Cluster'][concat_fig_list.index(fig)] = cluster_ID

    # Get pubseed_fig_information for fig in fig_dict['Fig']
    fig_information_list = get_multiple_pubseed_fig_information(fig_dict['Fig'], workers)
    for fig_info in fig_information_list:
        for var in fig_info:
            try:
                fig_dict[var].append(fig_info[var])
            except KeyError:
                fig_dict[var] = [fig_info[var]]

    fig_df = pd.DataFrame(fig_dict)
    updated_concat_df = pd.merge(concat_df, fig_df, on='Fig')
    return updated_concat_df

def build_ref_profiling_vector(fig_object_dict,
                               genomes_fig_to_match):

    similarity_variables = [fig for fig_list in genomes_fig_to_match.values() for fig in fig_list]
    operon_variables = ['ID', 'functional_role', 'length', 'start', 'end', 'operon', 'cluster', 'peg', 'gap', 'ev', 'strand']
    similarity_var = 'average_normalized_bit_score'

    reference_dict = {}
    for var in similarity_variables:
        reference_dict[var] = []

    # Build reference_dict to turn into pd.DataFrame
    for fig in fig_object_dict:
        peg_candidates = fig_object_dict[fig].peg_candidates

        # add operon information
        fig_dict = fig_object_dict[fig].__dict__
        for var in operon_variables:
            if var not in reference_dict:
                reference_dict[var] = []
            reference_dict[var].append(fig_dict[var])

        # add similarity value
        for var in similarity_variables:
            if var in peg_candidates:
                value = peg_candidates[var][similarity_var] if var in peg_candidates else 0.0
                reference_dict[var].append(value)
            else:
                reference_dict[var].append(0.0)

    ref_df = pd.DataFrame(reference_dict)
    return ref_df, similarity_variables

def get_ref_profiling_vectors(ref_fig_object_dict,
                              Blast_parameters,
                              workers):

    print('get_ref_profiling_vectors')
    if ref_fig_object_dict and Blast_parameters and workers:
        ref_fig_object_dict, genomes_fig_to_match = build_seq_similarity_network(ref_fig_object_dict,
                                                                                 Blast_parameters,
                                                                                 workers)
        # fig_object_dict = pickle.load(open('Inputs/pickle/fig_object_dict.pkl', "rb"))
        genomes_fig_to_match = get_genomes_fig_to_match(ref_fig_object_dict)

        ref_df, variables = build_ref_profiling_vector(ref_fig_object_dict, genomes_fig_to_match)
        return ref_df, ref_fig_object_dict, variables, genomes_fig_to_match

def retreive_fig_self_similarity_for_genomes(ref_fig_object_dict,
                                             Blast_parameters,
                                             workers):

    fig_list = list(ref_fig_object_dict.keys())
    genome_ID_list = list(set([get_orgid_from_fig(fig) for fig in fig_list]))
    for genome_ID in genome_ID_list:
        ref_fig_object_dict = retreive_fig_self_similarity(genome_ID,
                                                           ref_fig_object_dict,
                                                           Blast_parameters,
                                                           workers)
    return ref_fig_object_dict

def get_reference_set(PubSeed_speadsheet_file,
                      gold_standard_proteins_file,
                      prot_list,
                      Blast_parameters,
                      cluster_threshold,
                      workers,
                      reaction_genes,
                      usearch_path
                      ):

    # Get Gold Standard Protein and PubSeed protein data
    pubseed_df = load_pubseed_spreadsheet_cluster_content(PubSeed_speadsheet_file,
                                                          prot_list
                                                          )
    gsp_df = gold_standard_proteins(gold_standard_proteins_file)

    # concatenate pubseed_df and gsp_df. Instantiate objects from this concatenate protein data
    updated_concat_df = retreive_pubseed_information(pubseed_df, gsp_df, cluster_threshold, workers)

    # pickle.dump(updated_concat_df, open('Inputs/pickle/updated_concat_df.pkl', "wb"),
    #             protocol=pickle.HIGHEST_PROTOCOL)
    # updated_concat_df = pickle.load(open('Inputs/pickle/updated_concat_df.pkl', "rb"))
    # print('updated_concat_df', len(updated_concat_df))

    ref_fig_object_dict = instantiate_reference_fig(updated_concat_df)
    # Reduce the reference dataset
    print('ref_fig_object_dict', len(ref_fig_object_dict))
    # pickle.dump(ref_fig_object_dict, open('Inputs/pickle/instantiate_reference_fig.pkl', "wb"),
    #             protocol=pickle.HIGHEST_PROTOCOL)
    # ref_fig_object_dict = pickle.load(open('Inputs/pickle/instantiate_reference_fig.pkl', "rb"))
    ref_fig_object_dict = reduce_refdata(ref_fig_object_dict,
                                         reaction_genes,
                                         usearch_path)
    print('ref_fig_object_dict', len(ref_fig_object_dict))

    # Updated fig information with its self similarity
    ref_fig_object_dict = retreive_fig_self_similarity_for_genomes(ref_fig_object_dict,
                                                                   Blast_parameters,
                                                                   workers)
    return ref_fig_object_dict

if __name__ == "__main__":
    tstart = time()
    print("Start")

    # Define parameters of the Blast
    cutoff = 1e-05
    Blast_parameters = {"seq_type":"aa",
              "filter" : "F",
              "evalue" : str(cutoff),
              "wsize" : "0",
              "fasta" : '',
              "organism" : '',
              "act" : "BLAST"
            }
    workers = 50
    cluster_threshold = 5

    if True:
        ''' load pubseed spreadsheet cluster content '''
        PubSeed_speadsheet_file = 'Inputs/PubSeed_speadsheets/Acetate production.tsv'
        prot_list = ['HdrA', 'HdrB', 'HdrC', 'MvhD']
        gold_standard_proteins_file = 'Inputs/Gold_standard_proteins/Ref_operons.txt'
        reaction_genes = 'heterodisulfide_reductase'
        usearch_path = os.path.join('..', '..', 'SW', 'Usearch', 'usearch10.0.240_win32.exe')

        ref_df, ref_fig_object_dict, variables, genomes_fig_to_match = get_reference_set(PubSeed_speadsheet_file,
                                                                                         gold_standard_proteins_file,
                                                                                         prot_list,
                                                                                         Blast_parameters,
                                                                                         cluster_threshold,
                                                                                         workers,
                                                                                         reaction_genes,
                                                                                         usearch_path
                                                                                         )

    tend = time()
    print('Process duration :',int(tend - tstart))
    print ("Done")