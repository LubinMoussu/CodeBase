# -*- coding: utf-8 -*-
"""
Created on Thu Sep  6 13:52:42 2018

@author: lubin.moussu
"""

from urllib.parse import urlencode
from urllib.request import Request, urlopen
from pprint import pprint
import re 
from os import startfile
from operator import itemgetter
from itertools import groupby
import os
from re import findall, match
import numpy as np, scipy.stats as st, pandas as pd
from collections import Counter
from html.parser import HTMLParser
from urllib.parse import urlparse, parse_qs
from random import seed, randint, SystemRandom, sample
from webbrowser import open as webbrowser_open
import pickle
from multiprocessing.dummy import Pool as ThreadPool
import itertools
import urllib
from subprocess import Popen
from time import sleep
import operator
from http.client import *

from selenium import webdriver
from selenium.webdriver.firefox.firefox_binary import FirefoxBinary
from selenium.webdriver.common.keys import Keys
import selenium.webdriver.support.ui as ui
from selenium.webdriver.firefox.options import Options

from similaritymeasures import jaccard_similarity
from selenium_browser import *
from Multi_project_scripts import get_pubseed_single_fig_information

import bs4
try: 
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

    
class Protein_Evidence():
    def __init__(self, protein_ID,protein_IDs,genome_ID,evidence,pubmed_IDs):
        self.protein_ID = protein_ID
        self.protein_IDs = protein_IDs
        self.genome_ID = genome_ID
        self.evidence = evidence
        self.pubmed_IDs = pubmed_IDs

class Protein():
    def __init__(self,
                 accession_code,
                 database_sources,
                 entry_ID,
                 functional_roles,
                 genome_names,
                 matches_PubSeed_features,
                 pubmed_IDs,
                 genes
                 ):

        self.accession_code = accession_code
        self.database_sources = database_sources
        self.entry_ID = entry_ID
        self.functional_roles = functional_roles
        self.genome_names = genome_names
        self.matches_PubSeed_features = matches_PubSeed_features
        self.pubmed_IDs = pubmed_IDs
        self.genes = genes

    def get_PubSeed_Matches_features_from_database_Sources(self, Blast_parameters, all_PubSeed_genomes):

        entries = self.accession_code
        prot = False
        self.matches_PubSeed_features = []
        for entry in entries:
            matches_features = get_PubSeed_Matches_features_from_ProteinID(entry)
            if not matches_features:
                # If entry_ID of databases are not found in PubSeed, the protein sequence is retrieved from databases and blast against genomes of PubSeed
                uniprot_entry = get_uniprot_fasta(entry)
                if uniprot_entry and uniprot_entry['taxon_identifier']:
                    # If there is a identified reference genome, find PubSeed candidates
                    matches_features = get_PubSeed_Matches_features_from_Sources_prot_sequences(entry, Blast_parameters, uniprot_entry['seq'], uniprot_entry['taxon_identifier'], all_PubSeed_genomes)
                    
            if matches_features:
                # Check if genome_ID is within expected genome_IDs in all_PubSeed_genomes
                matches_features = [fig for fig in matches_features if get_orgid_from_fig(fig) in all_PubSeed_genomes]
                if matches_features:
                    self.matches_PubSeed_features = matches_features
                    prot = True

                    # If there are matches, add the gene names from Uniprot
                    try:
                        uniprot_entry
                    except NameError:
                        uniprot_entry = get_uniprot_fasta(entry)
                    if uniprot_entry and uniprot_entry['genes']:
                        self.genes = uniprot_entry['genes']
        return self, prot

    def update_literature_pubseed_feature(self, browser, EC_number_sufficiency):

        entry_pubmed_IDs = list(set([pubmed_ID for site in self.pubmed_IDs for pubmed_ID in self.pubmed_IDs[site]]))
        protein_ID = self.entry_ID
        features_to_curate = []
        consistency = False

        if protein_ID and self.matches_PubSeed_features and entry_pubmed_IDs:
            # if the protein_ID has an entry in PubSeed and an evidence article, check if the article is assigned to PubSeed entry
            for fig in self.matches_PubSeed_features:
                all_fig_information = get_pubseed_single_fig_information(fig)
                # Curate literature of the fig
                curate_literature_to_fig(fig, all_fig_information['Ev'], entry_pubmed_IDs, browser)
                print('{0} pubmed_ID_list was added to {1}, {2}'.format(str(entry_pubmed_IDs),
                                                                        str(fig),
                                                                        str(protein_ID)
                                                                        ))

                # check if the current annotation is consistent with reference annotations. If not, trigger manual curation
                consistency = functional_role_consistency(all_fig_information, self, 0.6, EC_number_sufficiency)

                if not consistency:
                    features_to_curate.append(fig)
        return consistency, features_to_curate

def PaperBlast_Query(seq):

        url =  'http://papers.genomics.lbl.gov/cgi-bin/litSearch.cgi?'
        values = {"query": seq,
                  "Search": 'Search'
                }

        for key in values:
            assert type(values[key]=="str")
        
        # Encode values in url code
        data = urlencode(values)
        # Create the request
        request = "{0}&{1}".format(url, data)
        opener = urllib.request.build_opener()
        try:
            response = opener.open(request, timeout=60)
            sleep(5)
            print('Sucessful PaperBlast_Query')
            return response
        except Exception as ex:
            print('PaperBlast_Query', ex)
            print(seq)
            return None
    
def get_entry_evidences(seq, targeted_sites):
    '''
    '''
    assert type(seq)==str
    page = PaperBlast_Query(seq)
    if page:
        soup = BeautifulSoup(page, 'html.parser')

        evidence_sets = {}
        genomes_features = soup.find_all('p', attrs={'style':"margin-top: 1em; margin-bottom: 0em;"})
        for genome in genomes_features:
            genome_ID = genomes_features.index(genome)

            entry_evidence = {'entry_ID': '',
                              'database_sources': {},
                              'accession_code': [],
                              'pubmed_IDs': [],
                              'functional_roles': {},
                              'genome_names': [],
                              'genes': []
                              }
            Accession_Code = ''
            functional_role = ''
            pubmed_IDs = []
            database_name = ''

            database_sources = {site: [] for site in targeted_sites}
            databases_pubmed_IDs = {site: [] for site in targeted_sites}
            databases_functional_role = {site: [] for site in targeted_sites}

            if genome_ID:
                for child in genome.children:
                    if child.name:
                        if child.name == 'br' and Accession_Code and pubmed_IDs:
                            # add information to entry_evidence corresponding to each database if the entry is a reference protein
                            database_sources[database_name].append(Accession_Code)
                            databases_pubmed_IDs[database_name].extend(pubmed_IDs)
                            databases_functional_role[database_name].append(functional_role)

                            Accession_Code = ''
                            functional_role = ''
                            pubmed_IDs = []

                        # database_entries
                        if child.has_attr('title') and child['title'] in targeted_sites:
                            database_name = str(child['title'])
                            # print(child.contents, child['title'])
                            raw_accession_Code = child.contents[0]
                            Accession_Code = raw_accession_Code if '/' not in raw_accession_Code \
                                else raw_accession_Code.split(' / ')[1]
                            Accession_Code = str(Accession_Code)

                            if '/' in raw_accession_Code:
                                genes = raw_accession_Code.split('_')[0] if '_' in raw_accession_Code\
                                    else raw_accession_Code.split(' / ')[0]
                                entry_evidence['genes'].append(genes)

                        # functional_role
                        if child.name == 'b':
                            functional_role = str(child.contents[0])
                            # functional_role = str(child.contents[0]).encode('utf-8')

                        # Genome_name
                        if child.name == 'i':
                            genome_name = str(child.string)
                            entry_evidence['genome_names'].append(genome_name)

                        # Pubmed entries
                        if child.has_attr('onmousedown'):
                            pubmed_url = child['href']
                            if pubmed_url and Accession_Code:
                                pubmed_IDs = pubmed_url.lstrip('http://www.ncbi.nlm.nih.gov/pubmed/').split(',')
                                pubmed_IDs = [str(elt) for elt in pubmed_IDs if elt.isdigit()]

            protein_IDs = list(set([protein_ID for site in database_sources for protein_ID in database_sources[site] if
                                    protein_ID[0].isupper()]))

            if protein_IDs:
                # For all protein_ID founds in protein_IDs, an object is instantiated with all features collected above
                for protein_ID in protein_IDs:
                    if protein_ID not in evidence_sets:
                        if entry_evidence['genome_names']:
                            matches_PubSeed_features = []
                            evidence_sets[protein_ID] = Protein(protein_IDs,
                                                                database_sources,
                                                                protein_ID,
                                                                databases_functional_role,
                                                                entry_evidence['genome_names'],
                                                                matches_PubSeed_features,
                                                                databases_pubmed_IDs,
                                                                entry_evidence['genes']
                                                                )
        return evidence_sets
   
def query_protein_ID_to_PubSeed(protein_ID):
    url = "http://pubseed.theseed.org//seedviewer.cgi?pattern={0}&page=Find&act=check_search".format(protein_ID)
    # Check if the protein_ID coming from different databases is present in PubSeed
    request = url
    opener = urllib.request.build_opener()
    try:
        response = opener.open(request, timeout=30)
        sleep(1)
    except urllib.error.URLError:
        return None
    return response

def get_orgid_from_fig(fig):
    org_id = findall('\d+\.\d+', fig)[0]
    return org_id

def get_uniprot_fasta(uniprot_ID):
    url = "https://www.uniprot.org/uniprot/{0}.fasta".format(uniprot_ID)
    request = url
    opener = urllib.request.build_opener()
    try:
        page = opener.open(request, timeout=120)
        sleep(1)
    # except (urllib.error.URLError, http.client.RemoteDisconnected) as error:
    #     return None
    except Exception as ex:
        print('uniprot', uniprot_ID, ex)

    response = [elt.decode("utf-8").rstrip("\n") for elt in page if elt.decode("utf-8").rstrip("\n")]
    seq = ''.join(response[1:])

    if response:
        ID = response[0]

        OS = findall('OS=(.+) OX=', ID)[0]
        taxon_identifier = findall(' OX=(\d*) ', ID)[0]
        genes = findall(' GN=(\w+)', ID)

        if 'strain' in OS:
            taxo = OS.split(' (')
            genus = taxo[0]
            strains = taxo[1].rstrip(')').lstrip('strain ').split(' / ')
            genome_names = []
            for strain in strains:
                genome_name = '{0} {1}'.format(genus, strain)
                genome_names.append(genome_name)
        else:
            genome_names = [OS]

        uniprot_entry = {'ID': uniprot_ID,
                         'seq' :seq,
                         'genome_names': genome_names,
                         'taxon_identifier': taxon_identifier,
                         'genes': genes
                         }
        return uniprot_entry

def get_PubSeed_Matches_features_from_ProteinID(protein_ID):

    page = query_protein_ID_to_PubSeed(protein_ID)
    soup = BeautifulSoup(page, 'html.parser')

    matches_features = []
    for link in soup.find_all('input', attrs={'id': "table_data_0"}):
        figs = findall('fig\|\d+\.\d+\.peg.\d+', str(link))
        for fig in figs:
            matches_features.append(fig)
    matches_features = list(set(matches_features))
    return matches_features

def get_PubSeed_genomes_with_Taxon_identifier(all_PubSeed_genomes,
                                              Taxon_identifier):

    matches = [elt for elt in all_PubSeed_genomes if Taxon_identifier in elt]
    return matches

def get_PubSeed_Matches_features_from_Sources_prot_sequences(id_seq,
                                                             Blast_parameters,
                                                             fasta_seq,
                                                             Taxon_identifier,all_PubSeed_genomes):
    '''
    Aim: protein sequence is blasted against PubSeed genome to get the right protein
    :param id_seq:
    :param Blast_parameters:
    :param fasta_seq:
    :param genome_ID:
    :return:
    '''

    matches_genomes = get_PubSeed_genomes_with_Taxon_identifier(all_PubSeed_genomes,Taxon_identifier)
    candidate_figs = []
    for genome_ID in matches_genomes:
        id_seq, peg_evalue, genome_ID = Blast_query_PEG_evalues(id_seq, Blast_parameters, fasta_seq, genome_ID)

        for candidate_fig in peg_evalue:
            if peg_evalue[candidate_fig][0] == 0.0:
                candidate_figs.append(candidate_fig)
    return candidate_figs

def get_PubSeed_Matches_features_from_Sources_IDs(evidence_sets,
                                                  Blast_parameters,
                                                  all_PubSeed_genomes,
                                                  workers):
    '''
    Aim : Find the PubSeed Fig associated with the protein_ID. Only keep protein_IDs for which PubSeed entries have been found
    :param evidence_sets:
    :param Blast_parameters:
    :param all_PubSeed_genomes:
    :return: evidence_sets
    '''

    workers = 50
    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results

    print(len(evidence_sets))

    args = zip(list(evidence_sets.values()),
               itertools.repeat(Blast_parameters),
               itertools.repeat(all_PubSeed_genomes)
               )

    results = pool.starmap((lambda x, y, z: x.get_PubSeed_Matches_features_from_database_Sources(y,
                                                                                                 z)), args)
    # results = pool.starmap(update_refPEGs_from_candidate_fig, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    for protein in results:
        protein_object, prot = protein
        # If the protein_ID cannot be found in PubSeed, the protein_object is removed
        if not prot:
            try:
                del evidence_sets[protein_object.entry_ID]
            except KeyError:
                pass

    return evidence_sets

def EC_number_similarity(EC1, EC2):

    EC1 = EC1.split('.')
    EC2 = EC2.split('.')

    si = [1 if EC1[i]==EC2[i] else 0 for i in range(len(EC1))]

def uniprot_potential_functional_roles(reference_roles):

    all_uniprot_names = []
    for reference_role in reference_roles:
        names = reference_role.split('; ')
        EC_number_pattern = 'EC (\d+.\d+.[\d-]+.[\d-]+)'
        EC_number = [elt for elt in names if findall(EC_number_pattern,elt)]
        if EC_number:
            EC_number = EC_number[0]
            names = ['{0} ({1})'.format(elt, EC_number) for elt in names if elt!=EC_number and len(elt)>5]
            all_uniprot_names.extend(names)
    return list(set(all_uniprot_names))

def update_functional_role_consistency(current_role_EC,
                                       reference_role_EC,
                                       consistency,
                                       current_role,
                                       reference_role,
                                       jaccard_threshold,
                                       EC_number_sufficiency):

    if current_role_EC and reference_role_EC:
        reference_role_EC = reference_role_EC[0] if reference_role_EC else ''
        if reference_role_EC == current_role_EC:
            # If both EC numbers are equal,
            if EC_number_sufficiency:
                consistency = True
            else:
                jaccard_score = jaccard_similarity(current_role, reference_role)
                if jaccard_score >= jaccard_threshold:
                    consistency = True
    return consistency

def functional_role_consistency(fig_information,
                                evidence_set,
                                jaccard_threshold,
                                EC_number_sufficiency):
    '''
    Aim: check that the functional role fits ontology systems role2
    Process: (1) take role1, role2
             (2) if EC_number_sufficiency, just check that EC number is the same
             (3) if not EC_number_sufficiency, calculate jaccard_similarity between both roles
             (4) return bolean of role1 role2 consistency
    :param role1:
    :param role2:
    :return:
    '''
    EC_number_pattern = 'EC (\d+.\d+.[\d-]+.[\d-]+)'
    current_role = fig_information['Function']
    current_role_EC = findall(EC_number_pattern,current_role)
    current_role_EC = current_role_EC[0] if current_role_EC else ''
    reference_roles = evidence_set.functional_roles

    consistency = False
    for site in reference_roles:
        if not consistency:
            if site == 'SwissProt':
                list_reference_roles = uniprot_potential_functional_roles(reference_roles[site])
            else:
                list_reference_roles = list(reference_roles[site])

            for reference_role in list_reference_roles:
                if not consistency:
                    reference_role_EC = findall(EC_number_pattern, reference_role)
                    consistency = update_functional_role_consistency(current_role_EC,
                                                                       reference_role_EC,
                                                                       consistency,
                                                                       current_role,
                                                                       reference_role,
                                                                       jaccard_threshold,
                                                                       EC_number_sufficiency)
    return consistency

def curate_literature_to_fig(fig,
                             current_pubmed_IDs,
                             entry_pubmed_IDs,
                             browser):

    # if the PubSeed fig is not associated with any pubmed_ID which however should be associated, add the pubmed_ID to PubSeed_ID
    if not current_pubmed_IDs:
        # Define pubmed_ID to add to the fig
        pubmed_IDs_to_add = [pubmed_ID for pubmed_ID in entry_pubmed_IDs if pubmed_ID not in current_pubmed_IDs]

        for pubmed_ID in pubmed_IDs_to_add:
            add_pubmedID_to_fig_selenium(fig, pubmed_ID, browser)
            sleep(3)
            break

def update_literature_PubSeed(evidence_sets,
                              EC_number_sufficiency):

    '''
    Aim : Automatically update literature in PubSeed for figs from evidence_sets
    :param evidence_sets:
    :return: figs for which current_role is different from reference_role
    '''

    # init_headless_firefox_browser
    browser = init_headless_firefox_browser()

    # init_pubseed_login
    ID, password =  'moussulubin', 'E2UGSat8z<'
    browser = init_pubseed_login(browser, ID, password)
    sleep(3)
    print('Connected to PubSeed')

    protein_features_to_curate = {}
    for protein_ID in list(evidence_sets.keys()):
        consistency, features_to_curate = evidence_sets[protein_ID].update_literature_pubseed_feature(browser, EC_number_sufficiency)
        if not consistency and features_to_curate:
            protein_features_to_curate[protein_ID] = list(set(features_to_curate))
    browser.quit()
    return protein_features_to_curate

def get_multiple_entry_evidences(ref_fig_list,
                                 targeted_sites):

    multiple_evidence_sets = {}
    fig_list = []
    for elt in ref_fig_list:
        seq = elt['Seq']
        fig_list.append(elt['Fig'])
        evidence_sets = get_entry_evidences(seq, targeted_sites)
        if evidence_sets:
            multiple_evidence_sets.update(evidence_sets)
    return multiple_evidence_sets, fig_list

def get_multithreading_multiple_entry_evidences(ref_fig_list,
                                                targeted_sites,
                                                workers):

    # make the Pool of workers
    pool = ThreadPool(workers)
    print('workers', workers)
    # and return the results

    fig_list = []
    seq_list = []
    for elt in ref_fig_list:
        seq_list.append(elt['Seq'])
        fig_list.append(elt['Fig'])

    args = zip(seq_list,
               itertools.repeat(targeted_sites))

    results = pool.starmap(get_entry_evidences, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    multiple_evidence_sets = {}
    for evidence_sets in results:
        if evidence_sets:
            multiple_evidence_sets.update(evidence_sets)
    return multiple_evidence_sets, fig_list

def get_arborescence_PubSeed_Matches_features_from_Sources_IDs(ref_fig_list,
                                                               matches_per_entry,
                                                               targeted_sites,
                                                               Blast_parameters,
                                                               all_PubSeed_genomes,
                                                               depth,
                                                               workers):

    iter = 0
    multiple_arborescence_evidence_sets = {}
    already_explored_fig_list = []

    pprint(matches_per_entry)
    for protein_ID in matches_per_entry:
        protein_IDs = database_sources = genome_names = genes = []
        databases_functional_role = databases_pubmed_IDs = {}

        multiple_arborescence_evidence_sets[protein_ID] = Protein(protein_IDs,
                                                                  database_sources,
                                                                  protein_ID,
                                                                  databases_functional_role,
                                                                  genome_names,
                                                                  [matches_per_entry[protein_ID]],
                                                                  databases_pubmed_IDs,
                                                                  genes)
        already_explored_fig_list.append(matches_per_entry[protein_ID])

    while iter < depth:
        evidence_sets, fig_list = get_multithreading_multiple_entry_evidences(ref_fig_list, targeted_sites, workers)
        already_explored_fig_list.extend(fig_list)

        # get PubSeed Matches only for protein_ID not present in multiple_arborescence_evidence_sets
        evidence_sets = {key:value for key,value in evidence_sets.items() if key not in multiple_arborescence_evidence_sets}

        evidence_sets = get_PubSeed_Matches_features_from_Sources_IDs(evidence_sets,
                                                                      Blast_parameters,
                                                                      all_PubSeed_genomes,
                                                                      workers)

        print('{0} reference external links associated with PubSeed features'.format(len(evidence_sets.keys())))

        # update multiple_arborescence_evidence_sets
        multiple_arborescence_evidence_sets.update(evidence_sets)

        # define new fig to find papers about it and its homologs
        new_ref_fig_list = [fig for protein_ID in evidence_sets for fig in evidence_sets[protein_ID].matches_PubSeed_features]
        new_ref_fig_list = [fig for fig in new_ref_fig_list if fig not in ref_fig_list and fig not in already_explored_fig_list]
        ref_fig_list = get_multiple_fasta_information(new_ref_fig_list)

        iter+=1
    return multiple_arborescence_evidence_sets

def write_evidence_sets_PubSeed_Matches(evidence_sets,
                                        gsp_file):

    gsp_dict = {'Peg': [],
                'Fig': []}

    for protein_ID in evidence_sets:

        peg = evidence_sets[protein_ID].genes[0] if evidence_sets[protein_ID].genes else protein_ID
        for fig in evidence_sets[protein_ID].matches_PubSeed_features:
            peg = peg.title()
            gsp_dict['Peg'].append(peg)
            gsp_dict['Fig'].append(fig)

    ref_df = pd.DataFrame(gsp_dict)
    ref_df.set_index('Peg', inplace=True)

    ref_df.to_csv(gsp_file, sep='\t', encoding='utf-8')
    return ref_df

def get_PubSeed_Matches_features_from_protein_IDs(reaction_genes,
                                                  protein_IDs_to_get,
                                                  targeted_sites,
                                                  Blast_parameters,
                                                  all_PubSeed_genomes,
                                                  depth,
                                                  workers,
                                                  EC_number_sufficiency
                                                  ):

    # get protein_IDs that do not match the PubSeed_ID format to search them in PubSeed
    protein_IDs_in_PubSeed = [elt for elt in protein_IDs_to_get
                              if findall("(fig.\d *.\d *.peg.\d *)", elt)
                              ]
    protein_IDs_not_fig_list = [elt for elt in protein_IDs_to_get
                                if elt not in protein_IDs_in_PubSeed
                                ]
    fig_list = []
    if protein_IDs_not_fig_list:
        # Query protein_IDs on PubSeed and get the most recently updated genome [-1]
        matches_per_entry = {entry: get_PubSeed_Matches_features_from_ProteinID(entry)
                             for entry in protein_IDs_to_get}
        matches_per_entry = {entry: matches[-1]
                             for entry, matches in matches_per_entry.items()
                             if matches}
        fig_list = [values for values in matches_per_entry.values()]
    fig_list.extend(protein_IDs_in_PubSeed)

    # Get evidence_sets for all reference_fig in fig_list
    print('{0} PubSeed entries found'.format(len(fig_list)))
    ref_fig_list = get_multiple_fasta_information(fig_list)
    evidence_sets = get_arborescence_PubSeed_Matches_features_from_Sources_IDs(ref_fig_list,
                                                                               matches_per_entry,
                                                                               targeted_sites,
                                                                               Blast_parameters,
                                                                               all_PubSeed_genomes,
                                                                               depth,
                                                                               workers)

    gsp_file = 'Inputs/Gold_standard_proteins/{0}.txt'.format(reaction_genes)
    # file_out = 'Inputs/PubSeed_matches/{0}.txt'.format(reaction_genes)
    ref_df = write_evidence_sets_PubSeed_Matches(evidence_sets, gsp_file)

    protein_features_to_curate = update_literature_PubSeed(evidence_sets, EC_number_sufficiency)
    # protein_features_to_curate contains PubSeed fig that should be manually (1) checked (2) annotated
    return evidence_sets, protein_features_to_curate

if __name__ == "__main__":
    import time
    from Multi_project_scripts import *
    start = time.time()
    print("Start")

    ''' Parse args '''
    import argparse
    parser = argparse.ArgumentParser(description="build reference_set from PubSeed and GSPs")
    parser.add_argument('-l', '--list', help='str(list of proteins of interest)', type=str, required=True)
    parser.add_argument("-d", "--depth", help='depth defines the number of times newly found Gold Standard Proteins are re-used to find new homologs', type=int, required=True)
    parser.add_argument("-p", "--pathway", help='name of pathway', type=str)
    args = parser.parse_args()
    protein_IDs_to_get = [str(item) for item in args.list.split(',')]
    depth = args.depth
    reaction_genes = args.pathway

    from sys import getrecursionlimit, setrecursionlimit
    getrecursionlimit()
    max_rec = 0x10000000
    setrecursionlimit(max_rec)

    # targeted_sites = ['SwissProt', 'EcoCyc', 'MetaCyc', 'BRENDA', 'MicrobesOnline']
    targeted_sites = ['SwissProt', 'EcoCyc', 'MetaCyc', 'BRENDA']

    all_PubSeed_genomes = get_csv_reader('Inputs/All_PubSeed_organisms.txt', delimiter=";")[0] # Only contains Bacteria + Archaea
    print('all_PubSeed_genomes', len(all_PubSeed_genomes))

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

    depth = 2  # Gold Standard Proteins found with 1 iteration can be used as inputs to find their homologs.
    # depth defines the number of times newly found Gold Standard Proteins are re-used to find new homologs
    workers = 3
    EC_number_sufficiency = True

    ''' Define proteins for which we want to find Gold Standard Proteins and their homologs
    '''
    protein_IDs_to_get = ['EG10585', 'Moth_1192', 'Moth_1193', 'Moth_1194', 'Moth_1195', 'Moth_1196']
    reaction_genes = 'heterodisulfide_reductase'  # name of the pathway for the storage of processing files

    if protein_IDs_to_get:
        evidence_sets, protein_features_to_curate = get_PubSeed_Matches_features_from_protein_IDs(reaction_genes,
                                                                                                  protein_IDs_to_get,
                                                                                                  targeted_sites,
                                                                                                  Blast_parameters,
                                                                                                  all_PubSeed_genomes,
                                                                                                  depth,
                                                                                                  workers,
                                                                                                  EC_number_sufficiency
                                                                                                  )

        pickle.dump(evidence_sets,
                    open('Inputs/pickle/{0}_evidence_sets.pkl'.format(reaction_genes), "wb"),
                    protocol=pickle.HIGHEST_PROTOCOL)
        print('{0} dump to {1}'.format('evidence_sets', 'Inputs/pickle/{0}_evidence_sets.pkl'.format(reaction_genes)))

        pickle.dump(protein_features_to_curate,
                    open('Inputs/pickle/protein_features_to_curate_{0}.pkl'.format(reaction_genes), "wb"),
                    protocol=pickle.HIGHEST_PROTOCOL)

    else:
        evidence_sets = pickle.load(
            open('Inputs/pickle/{0}_evidence_sets.pkl'.format(reaction_genes), "rb"))
        protein_features_to_curate = pickle.load(
            open('Inputs/pickle/protein_features_to_curate_{0}.pkl'.format(reaction_genes), "rb"))

        # pprint(evidence_sets)
        pprint(protein_features_to_curate)
        # TODO(take care of protein_features_to_curate)

    end = time.time()
    print(end - start)
    print ("Done")
