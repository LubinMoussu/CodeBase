
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
from time import sleep, time
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from webbrowser import open as webbrowser_open

from pandas_ml import ConfusionMatrix
import numpy as np
import pandas as pd
import scipy.stats as st

# from sklearn.naive_bayes import GaussianNB
# from sklearn.naive_bayes import ComplementNB
# from sklearn.model_selection import train_test_split
# from sklearn.metrics import confusion_matrix

from Multi_project_scripts import get_csv_reader, \
    get_orgid_from_fig, \
    merge_two_dicts, \
    is_orthologuous
from identifycandidates import identify_candidates_in_genomes
from evaluatecandidates import evaluate_candidate_cluster_list
from buildref import get_reference_set

def identify_evaluate_candidate_clusters_on_genome(genome_ID,
                                                   ref_fig_object_dict,
                                                   prot_list,
                                                   Blast_parameters,
                                                   workers,
                                                   thresholds):

    # Identify candidate cluster that could have proteins with same functions as proteins in prot_list
    # Homologs of proteins in prot_list are also searched in the candidate genome for further analysis

    ref_fig_object_dict, \
    candidate_fig_object_dict = identify_candidates_in_genomes(genome_ID,
                                                               ref_fig_object_dict,
                                                               prot_list,
                                                               Blast_parameters,
                                                               workers,
                                                               thresholds)

    # ref_fig_object_dict = pickle.load(open('Inputs/pickle/ref_fig_object_dict.pkl', "rb"))
    # candidate_fig_object_dict = pickle.load(open('Inputs/pickle/candidate_fig_object_dict.pkl', "rb"))

    all_candidate_annot = {}
    # candidate_fig_object_dict = pickle.load(open('Inputs/pickle/candidate_fig_object_dict.pkl', "rb"))
    if candidate_fig_object_dict:

        # get the determined functional roles of each candidate_fig
        candidate_fig_object_dict, \
        all_candidate_annot = evaluate_candidate_cluster_list(candidate_fig_object_dict,
                                                              ref_fig_object_dict,
                                                              thresholds)
    return genome_ID, candidate_fig_object_dict, all_candidate_annot

def identify_evaluate_candidate_clusters_on_multiple_genomes(list_organisms,
                                                             ref_fig_object_dict,
                                                             prot_list,
                                                             Blast_parameters,
                                                             workers,
                                                             thresholds):

    predictions = map(identify_evaluate_candidate_clusters_on_genome,
                      list_organisms,
                      itertools.repeat(ref_fig_object_dict),
                      itertools.repeat(prot_list),
                      itertools.repeat(Blast_parameters),
                      itertools.repeat(workers),
                      itertools.repeat(thresholds)
                      )
    predictions = [elt for elt in predictions]
    return predictions

def annotate_multiple_genomes(PubSeed_speadsheet_file,
                                            gold_standard_proteins_file,
                                            prot_list,
                                            Blast_parameters,
                                            thresholds,
                                            workers,
                                            reaction_genes,
                                            list_organisms):

    # Build reference set
    ref_fig_object_dict = get_reference_set(PubSeed_speadsheet_file,
                                            gold_standard_proteins_file,
                                            prot_list,
                                            Blast_parameters,
                                            thresholds['cluster_threshold'],
                                            workers,
                                            reaction_genes
                                            )

    # pickle.dump(ref_fig_object_dict, open('Inputs/pickle/ref_fig_object_dict.pkl', "wb"),
    #             protocol=pickle.HIGHEST_PROTOCOL)
    # ref_fig_object_dict = pickle.load(open('Inputs/pickle/ref_fig_network_object_dict.pkl', "rb"))

    if ref_fig_object_dict:
        for fig in ref_fig_object_dict:
            ref_fig_object_dict[fig].peg = ref_fig_object_dict[fig].peg.title()

        # Evaluate all the genomes on the list organisms
        predictions = identify_evaluate_candidate_clusters_on_multiple_genomes(list_organisms,
                                                                               ref_fig_object_dict,
                                                                               prot_list,
                                                                               Blast_parameters,
                                                                               workers,
                                                                               thresholds)

        return predictions

if __name__ == "__main__":

    tstart = time()
    print("Start")

    ''' Define all parameters
    '''
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
    workers = 10

    thresholds = {'cutoff': 1e-05,
                  'filter': 1e-10,
                  'operon_reannotation': 1e-35,
                  'prot_reannotation': 1e-275,
                  'prot_withdrawal': 1e-15,
                  'cluster_threshold': 5,
                  'orthologuous_thresholds': {'single': {'aa_identity': 0.4, 'length': 0.7},
                                             'cluster': {'aa_identity': 0.3, 'length': 0.6}
                                              }
        }

    Blast_parameters['cutoff'] = thresholds['cutoff']
    genome_ID = '339860.6'

    ''' Build ref
    '''
    PubSeed_speadsheet_file = 'Inputs/PubSeed_speadsheets/Acetate production.tsv'
    prot_list = ['HdrC', 'MvhD', 'HdrB', 'HdrA', 'MetF', 'metV']
    prot_list = [elt.title() for elt in prot_list]

    gold_standard_proteins_file = 'Inputs/Gold_standard_proteins/heterodisulfide_reductase.txt'
    reaction_genes = 'heterodisulfide_reductase'

    if prot_list:
        # List of all PubSeed genomes curated in Agora
        list_agora_organisms = get_csv_reader('Inputs/all_organisms.txt')[0][0].split(', ')

        predictions = annotate_multiple_genomes(PubSeed_speadsheet_file,
                                                gold_standard_proteins_file,
                                                prot_list,
                                                Blast_parameters,
                                                thresholds,
                                                workers,
                                                reaction_genes,
                                                list_agora_organisms)

            tend = time()
            print('Process duration : ', int(tend - tstart))
            print("Done")


