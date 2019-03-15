
from csv import reader as csv_reader
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from webbrowser import open as webbrowser_open
from pprint import pprint
from re import findall, match
from urllib.parse import urlparse, parse_qs
from time import sleep
from multiprocessing.dummy import Pool as ThreadPool
import itertools
from urllib.parse import urlparse, parse_qs

import webbrowser
import pickle
import urllib
import time
import bs4

try:
    from BeautifulSoup import BeautifulSoup
except ImportError:
    from bs4 import BeautifulSoup

from decorators import function_timing

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z

def Blast_query(Blast_parameters, genome_ID):

    '''
    Aim : Blast the seq into the genome_ID on PubSeed
    :param Blast_parameters:
    :param genome_ID:
    :return:
    '''

    url = "http://pubseed.theseed.org//seedviewer.cgi?page=BlastRun&"

    Blast_parameters['organism'] = genome_ID
    for key in Blast_parameters:
        assert type(Blast_parameters[key] == "str")

    data = urlencode(Blast_parameters)
    # Create the request
    request = "{0}&{1}".format(url, data)
    opener = urllib.request.build_opener()
    try:
        response = opener.open(request,timeout=60)
        sleep(1)
        return response
    except urllib.error.URLError or socket.timeout:
        return None

def read_fasta(file):
    with open(file, 'r') as fp:
        name, seq = None, []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if name: yield (name, ''.join(seq))
                name, seq = line, []
            else:
                line = line.replace(' ', '')
                line = line.rstrip('\n\r')
                seq.append(line)
        if name: yield (name, ''.join(seq))

def fasta_reader(file):
    ids_sequences = {}  # ids_sequences[id_seq] = {'Seqs':[list of sequences], 'Name':[names of protein], 'Fig':[]}
    sequence, ids, prot_name, fig = '', '', '', ''
    with open(file, "r") as f1:
        for line in f1:
            if line.startswith(">"):
                tot = line.split(' [')
                tot = tot[0] if len(tot) == 2 else ' ['.join(tot[:2])
                full_name = tot.lstrip('{0} '.format(ids)).rstrip(' ')

                if ids:
                    # add key to ids_sequences
                    if ids not in ids_sequences:
                        ids_sequences[ids] = {'Seqs': [], 'Names': [], 'Figs': []}

                    ids_sequences[ids]['Seqs'].append(sequence)
                    ids_sequences[ids]['Names'].append(prot_name)
                    ids_sequences[ids]['Figs'].append(fig)

                    ids = full_name.split(' ')[0].lstrip('>')
                    fig = findall('fig\|\d+\.\d+\.peg.\d+', full_name)
                    if fig:
                        fig = fig[0]
                        prot_name = full_name.split(fig + ' ')[1]
                    sequence = ""
                else:
                    ids = full_name.split(' ')[0].lstrip('>')
                    fig = findall('fig\|\d+\.\d+\.peg.\d+', full_name)
                    if fig:
                        fig = fig[0]
                        prot_name = full_name.split(fig + ' ')[1]
            else:
                sequence += line.rstrip("\n")

    if ids not in ids_sequences:
        ids_sequences[ids] = {'Seqs': [], 'Names': [], 'Figs': []}
    ids_sequences[ids]['Seqs'].append(sequence)
    ids_sequences[ids]['Names'].append(prot_name)
    if fig:
        ids_sequences[ids]['Figs'].append(fig)
    return ids_sequences

def get_csv_reader(file, delimiter="\t"):
    data=[]
    with open(file) as fd:
        for line in csv_reader(fd, delimiter=delimiter):
            data.append(line)
    return data

def Blast_PEG_evalues(response):
        ''' From a Blast file, get a dictionnary of PEGs and their related evalue
        '''
        peg_evalue = {}  # peg_evalue[peg] = {"length":, "Evalue":}
        try:
            response = [elt.decode("utf-8") for elt in response]

            for elt in response:
                if '<a href="?page=Annotation&feature=fig|' in elt and '><a name =' not in elt:
                    # Find the PEG id
                    expression = '">(.*?)</a>'
                    PEG = findall(expression, elt)[0]
                    elt = elt.replace(' ', '')

                    # Score =
                    expression = '(\d+>\d+)'
                    score = int(findall(expression, elt)[0].split('>')[1])

                    # Find the Evalue
                    subtable = elt.split('</a>')
                    assert len(subtable) == 3
                    Evalue = subtable[2]

                    if Evalue[0] == 'e':
                        Evalue = '1' + Evalue
                    Evalue = float(Evalue)
                    assert type(Evalue) == float

                    # save into the peg_evalue
                    if PEG and (Evalue or Evalue == 0.0):
                        peg_evalue[PEG] = [Evalue, score]
        except Exception as ex:
            print(ex)
        return peg_evalue

def get_PEG_evalues_from_Blast(response):

    page = response
    soup = BeautifulSoup(page, 'html.parser')
    matches_features = []

    for pre in soup.find_all('pre'):
        print(pre)

def Blast_query_PEG_evalues(id_seq, Blast_parameters, fasta_seq, genome_ID):
    '''
    :param id_seq:
    :param Blast_parameters:
    :param fasta_seq:
    :param genome_ID:
    :return:
    '''
    Blast_parameters['fasta'] = fasta_seq
    response = Blast_query(Blast_parameters, genome_ID)
    peg_evalue = Blast_PEG_evalues(response)

    return id_seq, peg_evalue, genome_ID

def get_fig_clusters(fig_list, cluster_threshold):

    positions_fig = {int(fig.split('.peg.')[1]): fig for fig in fig_list}
    positions = list(positions_fig.keys())
    cluster_groups = dict(enumerate(continuous_figs(positions, cluster_threshold), 1))
    for cluster in cluster_groups:
        figs_in_cluster = [positions_fig[elt] for elt in cluster_groups[cluster]]
        yield figs_in_cluster

def get_candidate_clusters_for_operon(peg_candidates, thresholds):

    candidate_figs = [k1 for k,v in peg_candidates.items() for k1 in v.keys()]
    return get_fig_clusters(candidate_figs, thresholds['cluster_threshold'])

def get_protein_significant_similarities(query_fig,
                                         Blast_parameters,
                                         fasta_seq,
                                         genome_ID):
    '''
    :param id_seq:
    :param Blast_parameters:
    :param fasta_seq:
    :param genome_ID:
    :return:
    '''
    Blast_parameters['fasta'] = fasta_seq
    page = Blast_query(Blast_parameters, genome_ID)

    soup = BeautifulSoup(page, 'html.parser')
    fig_similarity_information = {}
    if soup:
        for pre in list(soup.find_all('pre')):
            pre_children = [elt for elt in pre.children if elt.name == 'pre']
            for elt in pre_children[1:]:
                    fig_information = {}
                    children = [child for child in elt.children]
                    links = [child.name for child in children if child.name == 'a']
                    fig = ''
                    if len(links) == 2:
                        for child in children:
                            if child.name == 'a' and 'href' in child.attrs:
                                parsed_qs = parse_qs(child['href'])
                                if parsed_qs and 'feature' in parsed_qs.keys():
                                    fig = parsed_qs['feature'][0]
                            elif not child.name and fig:
                                # Check that child.string contains something less then space
                                if not all(element == ' ' for element in child.string):
                                    try:
                                        length = findall('Length\s*=\s*(\d*)', child.string)[0]
                                        score = findall('Score\s*=\s*\d*.\d*\s*bits\s*\((\d*)\)', child.string)[0]
                                        evalue = findall('Expect\s*=\s*(e-\d*|\de-\d*|\d.\d*),\s*Method', child.string)[0]
                                        evalue = '1{0}'.format(evalue) if evalue[0] == 'e' else evalue
                                        identities = findall('Identities\s*=\s*(\d*/\d*)\s*\((\d*)%\)', child.string)

                                        fig_information['length'] = int(length)
                                        fig_information['score'] = float(score)
                                        fig_information['evalue'] = float(evalue)
                                        try:
                                            fig_information['identities'] = identities[0]
                                        except IndexError:
                                            fig_information['identities'] = identities

                                    except Exception as ex:
                                        print('get_protein_significant_similarities', ex)
                                        print(query_fig, fig)
                                        print(type(child.string), len(child.string))
                                        print(links)
                                        print('---------------------------------------------')

                    if fig_information and fig:
                        fig_similarity_information[fig] = fig_information
    return query_fig, fig_similarity_information

def get_genomes_fig_to_match(fig_object_dict):

    '''
    Aim : For each genome, find reference proteins to match
    :param fig_object_dict:
    :return:
    '''
    genome_fig_dict = {}
    for fig in fig_object_dict:
        genome_ID = get_orgid_from_fig(fig)
        if genome_ID not in genome_fig_dict:
            genome_fig_dict[genome_ID] = []
        genome_fig_dict[genome_ID].append(fig)
    return genome_fig_dict

def get_protein_similarities_from_multithread(results):

    query_candidates = {}  # query_candidates[query][candidate] = similarities
    # Add subset to the peg_candidates

    for group in results:
        query, candidate_similarities = group
        if candidate_similarities:
            if query not in query_candidates:
                query_candidates[query] = candidate_similarities
            else:
                query_candidates[query].update(candidate_similarities)
    return query_candidates

def get_sequence_similarity_genome_with_refPEGs(genome_ID,
                                                figs_to_match,
                                                IDS_sequences,
                                                Blast_parameters,
                                                workers):

    assert type(figs_to_match) == list
    ''' Only keep similarity between reference proteins '''
    refFigs_candidates = get_multithreading_protein_significant_similarities(genome_ID,
                                                                             IDS_sequences,
                                                                             Blast_parameters,
                                                                             workers)

    def trim_refFigs_candidates(refFigs_candidates,
                                figs_to_match):
        # Del candidateFig that are not corresponding to reference Figs
        for candidate in list(refFigs_candidates.keys()):
            for candidateFig in list(refFigs_candidates[candidate].keys()):
                if candidateFig not in figs_to_match:
                    refFigs_candidates[candidate].pop(candidateFig)
        return refFigs_candidates

    refFigs_candidates = trim_refFigs_candidates(refFigs_candidates,
                                                 figs_to_match)
    return refFigs_candidates

def update_refPEGs(fig_object_dict,
                   reffigs_candidates):
    '''
    Aim: Update the fig_object_dict with the similarity network matching the genome_ID
    :param fig_object_dict:
    :param reffigs_candidates:
    :return:
    '''
    fig_candidate_list = list(set([fig for fig_list in reffigs_candidates.values() for fig in fig_list]))
    for query in reffigs_candidates.keys():
        for candidate_fig in reffigs_candidates[query]:
            if candidate_fig in fig_object_dict:
                fig_object = fig_object_dict[candidate_fig]
                fig_object.peg_candidates[query] = reffigs_candidates[query][candidate_fig]

    # update_average_normalized_bit_score to all fig_object
    fig_object_dict = update_average_normalized_bit_score(fig_object_dict)
    return fig_object_dict

@function_timing
def retreive_fig_self_similarity(genome_ID,
                                 fig_object_dict,
                                 Blast_parameters,
                                 workers):
    '''
    Aim: Blast the fig into its own genome to get the fig own similarity
    :param genome_ID:
    :param fig_object_dict:
    :param Blast_parameters:
    :param workers:
    :return:
    '''

    figs_sequences = {fig: [item.prot_sequence]
                      for fig, item in fig_object_dict.items()
                      }

    figs_to_match = list(fig_object_dict.keys())
    reffigs_candidates = get_sequence_similarity_genome_with_refPEGs(genome_ID,
                                                                     figs_to_match,
                                                                     figs_sequences,
                                                                     Blast_parameters,
                                                                     workers)
    fig_object_dict = update_refPEGs(fig_object_dict,
                                     reffigs_candidates)
    return fig_object_dict

def get_multithreading_protein_significant_similarities(genome_ID, ids_sequences, Blast_parameters, workers):

    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results
    id_seqs = []
    fasta_seqs = []
    for key, value in ids_sequences.items():
        for sub_value in value:
            id_seqs.append(key)
            fasta_seqs.append(sub_value)

    args = zip(id_seqs,
               itertools.repeat(Blast_parameters),
               fasta_seqs,
               itertools.repeat(genome_ID))

    results = pool.starmap(get_protein_significant_similarities, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    query_candidates = get_protein_similarities_from_multithread(results)
    return query_candidates

def get_query_fasta_sequences(feature):
    ''' Query fasta sequence of feature in PubSEED
    '''

    url = "http://pubseed.theseed.org//seedviewer.cgi?"
    values = {"page": "ShowSeqs",
              "feature": feature,
              "FLANKING": "500",
              "Sequence": "Protein Sequence",
              "Download": "Download Sequences",
              }
    for key in values:
        assert type(values[key] == "str")

    # Encode values in url code
    data = urlencode(values)
    # Create the request
    request = url + data
    opener = urllib.request.build_opener()

    while True:
        try:
            response = opener.open(request, timeout=60)
            sleep(1)
            return response
        except Exception as ex:
            print('query_fasta_sequences', feature, ex)

def query_fasta_sequences(feature):

    response = get_query_fasta_sequences(feature)
    response = [elt.decode("utf-8").rstrip("\n") for elt in response if elt.decode("utf-8").rstrip("\n")]
    seq = "".join(response[1:])
    full_name = response[0][1:]

    tot = full_name.split(' [')
    tot = tot[0] if len(tot) == 2 else ' ['.join(tot[:2])

    prot_name = tot.lstrip('{0} '.format(feature)).split('##')[0].rstrip(' ')
    feature_info = {"Fig": feature, "Prot_name": prot_name, "Seq": seq}
    return feature_info

def add_pubmedID_to_fig_selenium(fig, pubmed, browser):

    assert type(pubmed) == str
    assert pubmed.isdigit()
    url = 'http://pubseed.theseed.org/seedviewer.cgi?'

    values = {'page': 'EditPubMedIdsDSLits',
              'feature': fig,
              'PMID': pubmed,
              'SaveNew': 'Save to Dlits',
              'getPubmeds': '-1',
              'user': 'moussulubin'
              }

    # Encode values in url code
    fig = fig.encode(encoding='UTF-8', errors='strict')

    data = urllib.parse.urlencode(values)
    # Create the request
    request = "{0}{1}".format(url, data)

    try:
        browser.set_page_load_timeout(30)
        browser.get(request)
        print('{0} added to {1}'.format(pubmed, fig))
        return browser
    except Exception as ex:
        print('add_pubmedID_to_fig_selenium', ex)

def get_average_normalized_bit_score(fig_object1, fig_object2):

    peg_candidates1 = fig_object1.peg_candidates
    peg_candidates2 = fig_object2.peg_candidates

    fig1 = fig_object1.ID
    fig2 = fig_object2.ID

    similarity_score = peg_candidates1[fig2]['score'] if fig2 in peg_candidates1 else peg_candidates2[fig1]['score']
    self_score1 = peg_candidates1[fig1]['score'] if fig1 in peg_candidates1 else 0
    self_score2 = peg_candidates2[fig2]['score'] if fig2 in peg_candidates2 else 0

    if similarity_score and self_score1 and self_score2:
        average_normalized_bit_score = similarity_score / (float(self_score1 + self_score2)/2)
    else:
        average_normalized_bit_score = 0
    return average_normalized_bit_score

def get_fig_information(fig):
    url = 'http://pubseed.theseed.org//protein.cgi?'

    values = {'prot': fig,
              'user': 'moussulubin'
              }
    # Encode values in url code
    data = urlencode(values)
    # Create the request
    request = "{0}&{1}".format(url, data)
    opener = urllib.request.build_opener()
    while True:
        try:
            response = opener.open(request, timeout=120)
            sleep(1)
            return response
        except Exception as ex:
            print('get_fig_information', fig, ex)

def get_pubseed_fig_information(fig):

    page = get_fig_information(fig)
    soup = BeautifulSoup(page, 'html.parser')

    all_fig_information = {}
    for link in soup.body.find_all('span', attrs={'id':'chromosome_context_block_content'}):
        for table in link.find_all('table'):

            trs = list(table.find_all('tr'))
            header = [elt for elt in list(trs[0]) if elt.name == 'th']
            header = [elt.string if len(elt.contents)==1 else elt.contents[0] for elt in header]

            for tr in trs[1:]:
                feature_info = [elt for elt in list(tr) if elt.name == 'td']
                feature_info = [elt.contents if elt.find_all('a', href=True) else elt.string for elt in feature_info]

                fig_descripters = {var:feature_info[header.index(var)] for var in header}
                fig_descripters['Fid'] = parse_qs(urlparse(fig_descripters['Fid'][0]['href']).query)['feature'][0]

                if fig_descripters['Fid']:
                    fig_descripters['strand'] = fig_descripters.pop('\xa0')
                    # fig_descripters['Ev']
                    if fig_descripters['Ev'] and type(fig_descripters['Ev'][0]) == bs4.element.Tag:
                        fig_descripters['Ev'] = [str(a.string) for a in fig_descripters['Ev'][0].find_all('a', href=True)]
                    else:
                        fig_descripters['Ev'] = []

                    # fig_descripters['Function']
                    if fig_descripters['Function']:
                        if type(fig_descripters['Function']) == list:
                            fig_descripters['Function'] = "".join([elt.string for elt in fig_descripters['Function']])
                        else:
                            fig_descripters['Function'] = str(fig_descripters['Function'])
                    else:
                        fig_descripters['Function'] = ""

                    # clean
                    relevant_descripters = ['End','Ev','Fid','Function','Gap','Size','Start','strand']
                    fig_descripters = {key:fig_descripters[key] if type(fig_descripters[key]) == list else str(fig_descripters[key]) for key in relevant_descripters}

                    if fig_descripters['Fid'] not in all_fig_information:
                        all_fig_information[fig_descripters['Fid']] = fig_descripters
    return all_fig_information

def get_pubseed_single_fig_information(fig):
    '''
    Aim: From pubseed information from a whole cluster, get only the pubseed information for the fig
    :param fig:
    :return:
    '''
    return get_pubseed_fig_information(fig)[fig]

def get_multiple_pubseed_fig_information(fig_list, workers):

    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results

    args = zip(fig_list)
    results = pool.starmap(get_pubseed_single_fig_information(fig), args)

    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    return results

def continuous_figs(list_fig, threshold):
        assert type(list_fig) == list

        list_fig.sort()
        prev = None
        group = []
        for item in list_fig:
            if not prev or item - prev <= threshold:
                group.append(item)
            else:
                yield group
                group = [item]
            prev = item
        if group:
            yield group

def get_orgid_from_fig(fig):
    org_id = findall('\d+\.\d+', fig)[0]
    return org_id

def get_featureID_frim_fig(fig):

    featureID = findall('.peg.(\d*)', fig)[0]
    return featureID

def get_fig_from_genome_feature(genome_ID, feature_ID):

    fig = 'fig|{0}.peg.{1}'.format(genome_ID, feature_ID)
    return fig

def get_multiple_fasta_information(fig_list):

    for fig in fig_list:
        fig_info = query_fasta_sequences(fig)
        yield fig_info

def get_peg_candidates_from_multithread(results):
    peg_candidates = {}  # peg_candidates[peg][candidate] = evalue
    # Add subset to the peg_candidates

    for group in results:
        id_seq, peg_evalue, organism = group
        assert type(peg_evalue) == dict
        if peg_evalue and id_seq not in peg_candidates:
            peg_candidates[id_seq] = peg_evalue

        else:
            for k, v in peg_evalue.items():
                # If the fig is not already in peg_candidates, add it
                if k not in peg_candidates[id_seq]:
                    peg_candidates[id_seq][k] = v
                else:
                    # if there are already matches for this id_seq to this fig, replace e-value
                    current_Evalue = peg_candidates[id_seq][k][0]
                    new_Evalue = v[0]
                    if new_Evalue < current_Evalue:
                        # if current E.value < already saved E.value, change E.value
                        peg_candidates[id_seq][k][0] = new_Evalue
    return peg_candidates

def multithreaded(organism_id, ids_sequences, Blast_parameters, workers):
    # make the Pool of workers
    pool = ThreadPool(workers)
    # and return the results
    id_seqs = []
    fasta_seqs = []
    for key, value in ids_sequences.items():
        for sub_value in value:
            id_seqs.append(key)
            fasta_seqs.append(sub_value)

    args = zip(id_seqs,
               itertools.repeat(Blast_parameters),
               fasta_seqs,
               itertools.repeat(organism_id))

    results = pool.starmap(Blast_query_PEG_evalues, args)
    # close the pool and wait for the work to finish
    pool.close()
    pool.join()

    # Process multithreading process
    peg_candidates = get_peg_candidates_from_multithread(results)
    return peg_candidates

def is_orthologuous(fig_object1, fig_object2, orthologuous_thresholds):
    '''
    Aim : Check if fig1 is orthologous to fig2. Fig1 is the candidate and fig2 the query
    :param fig_object1:
    :param fig_object2:
    :param orthologuous_thresholds:
    :return:
    '''
    peg_candidates1 = fig_object1.peg_candidates
    peg_candidates2 = fig_object2.peg_candidates
    fig1 = fig_object1.ID
    fig2 = fig_object2.ID

    identities = ()
    query_length = 0

    if fig2 in peg_candidates1:
        identities = peg_candidates1[fig2]['identities']
        query_length = float(fig_object2.length)

    aa_identity = 0
    length_coverage = 0
    if identities and query_length:
        aa_identity = float(identities[1]) / 100
        alignment_length = float(identities[0].split('/')[1])
        length_coverage = round(alignment_length / query_length, 3)
        fig_object1.peg_candidates[fig2]['aa_identity'] = aa_identity
        fig_object1.peg_candidates[fig2]['length_coverage'] = length_coverage

        if aa_identity >= orthologuous_thresholds['aa_identity'] and length_coverage >= orthologuous_thresholds['length']:
            # update the fig object
            fig_object1.peg_candidates[fig2]['orthologuous'] = True
            return True, aa_identity, length_coverage
        fig_object1.peg_candidates[fig2]['orthologuous'] = False
    return False, aa_identity, length_coverage

def are_fig_contiguous(fig_list):

    cluster_list = list(get_fig_clusters(fig_list, 1))
    if len(cluster_list) == 1:
        return True
    return False

def query_best_hits(fig_object1, fig_object2, fig_object_dict):

    '''
    Aim : find best similarity hits in the targeted_genome from blast of the fig
    :param fig_object1:
    :return:
    '''
    genome1 = get_orgid_from_fig(fig_object1.ID)
    query = fig_object2.ID

    genome1_hits = [fig for fig in fig_object_dict if get_orgid_from_fig(fig) == genome1]
    genome1_bit_scores = list(set([fig_object_dict[fig1].peg_candidates[query]['average_normalized_bit_score']
                                             for fig1 in genome1_hits
                                             if query in fig_object_dict[fig1].peg_candidates.keys()
                                             and get_orgid_from_fig(query) != genome1
                                             ]))

    best_hits = []
    if genome1_bit_scores:
        genome1_bit_scores.sort(reverse=True)
        best_score = genome1_bit_scores[0]
        best_hits = [fig for fig in genome1_hits
                     if query in fig_object_dict[fig].peg_candidates
                     and fig_object_dict[fig].peg_candidates[query]['average_normalized_bit_score'] == best_score
                     ]
        return best_hits
    return best_hits

def is_best_hit(fig_object1, fig_object2, fig_object_dict):

    '''
    Aim : From the blast of fig2 against genome1, determine if fig1 is the best hit
    :param fig_object1:
    :param fig_object2:
    :return:
    '''

    fig1 = fig_object1.ID
    best_hits = query_best_hits(fig_object1, fig_object2, fig_object_dict)
    if fig1 in best_hits:
        return True
    return False

def are_bidirectional_best_hit(fig_object1, fig_object2, fig_object_dict):

    if is_best_hit(fig_object1, fig_object2, fig_object_dict) and \
            is_best_hit(fig_object2, fig_object1, fig_object_dict):
            return True
    return False

def are_TRIANGLE(fig_object1, fig_object2, fig_object3, fig_object_dict):
    '''
    Aim : xa in GA, xb in GB and xc in GC, form a TRIANGLE if and only if
    xa and xb are BBHs, xb and xc are BBHs and xc and xa are BBHs.
    :param fig_object1:
    :param fig_object2:
    :param fig_object3:
    :return:
    '''
    if are_bidirectional_best_hit(fig_object1, fig_object2, fig_object_dict) and \
        are_bidirectional_best_hit(fig_object1, fig_object3, fig_object_dict) and \
            are_bidirectional_best_hit(fig_object3, fig_object2, fig_object_dict):
        return True
    return False

def get_contiguous_orthologues(comb):
    # For each comb of the genome1, check if they are contiguous
    fig1_list = list(comb)
    if are_fig_contiguous(fig1_list):
        # if comb are contiguous, check that orthologues in genome2 are contiguous
        genome2_orthologues_match_comb = [fig2 for fig1 in fig1_list for fig2 in genome1_orthologues_dict[fig1]]
        if are_fig_contiguous(genome2_orthologues_match_comb):
            yield(fig1_list, genome2_orthologues_match_comb)

def gene_order_conservation(fig_object_dict1, fig_object_dict2, orthologuous_thresholds):
    '''
    Aim: Calculate the gene order conservation between 2 genomes reduced to a subsystem
    :param fig_object_dict:
    :return:
    process:
    1. get all fig from both genomes
    2. define orthologues
    3. get orthologues contiguous in each genome
    4. get couple of orthologues contiguous common on both
    '''

    fig_list_from_genome1 = [fig for fig in fig_object_dict1.keys()]
    fig_list_from_genome2 = [fig for fig in fig_object_dict2.keys()]

    combinations = list(itertools.product(fig_list_from_genome1, fig_list_from_genome2))
    orthologues = []
    genome1_orthologues_dict = {}

    cluster_comparison = True if (len(fig_list_from_genome1) > 1
                                  and len(fig_list_from_genome2) > 1) else False

    adapted_orthologuous_thresholds = orthologuous_thresholds['cluster'] \
        if cluster_comparison else orthologuous_thresholds['single']

    def get_genome1_orthologues_dict(combination):

        fig1 = combination[0]
        fig2 = combination[1]
        orthologuous, aa_identity, length_coverage = is_orthologuous(fig_object_dict1[fig1],
                                                                     fig_object_dict2[fig2],
                                                                     adapted_orthologuous_thresholds)
        if orthologuous:
            return combination

    orthologues = [elt for elt in map(get_genome1_orthologues_dict, combinations) if elt]
    genome1_orthologues_dict = {}
    for combination in orthologues:
        try:
            genome1_orthologues_dict[combination[0]].append(combination[1])
        except KeyError:
            genome1_orthologues_dict[combination[0]] = [combination[1]]


    genome1_orthologues = [comb[0] for comb in orthologues]
    genome1_ortho_combinations = list(itertools.combinations(genome1_orthologues, 2))

    # Find contiguous_orthologues in both clusters
    contiguous_orthologues = list(map(get_contiguous_orthologues, genome1_ortho_combinations))

    orthologue_number = len(orthologues)
    contiguous_orthologues_number = len(contiguous_orthologues)+1
    if orthologue_number:
        gene_order_conservation_score = float(contiguous_orthologues_number) / orthologue_number
        return gene_order_conservation_score, contiguous_orthologues_number, orthologue_number, orthologues
    return 0.0, contiguous_orthologues_number, orthologue_number, orthologues

def BBHs_order_conservation(fig_object_dict1, fig_object_dict2, fig_object_dict, orthologuous_thresholds):
    '''
    Aim: Calculate the gene order conservation between 2 genomes reduced to a subsystem
    :param fig_object_dict:
    :return:
    process:
    1. get all fig from both genomes
    2. define orthologues
    3. get orthologues contiguous in each genome
    4. get couple of orthologues contiguous common on both
    '''

    fig_list_from_genome1 = [fig for fig in fig_object_dict1.keys()]
    fig_list_from_genome2 = [fig for fig in fig_object_dict2.keys()]

    combinations = list(itertools.product(fig_list_from_genome1, fig_list_from_genome2))
    orthologues = []
    genome1_orthologues_dict = {}

    for combination in combinations:
        fig1 = combination[0]
        fig2 = combination[1]
        BBH = are_bidirectional_best_hit(fig_object_dict1[fig1], fig_object_dict2[fig2], fig_object_dict)
        if BBH:
            orthologues.append(combination)
            if combination[0] not in genome1_orthologues_dict:
                genome1_orthologues_dict[combination[0]] = []
            genome1_orthologues_dict[combination[0]].append(combination[1])

    genome1_orthologues = [comb[0] for comb in orthologues]
    genome1_ortho_combinations = list(itertools.combinations(genome1_orthologues, 2))

    # Find contiguous_orthologues in both clusters
    contiguous_orthologues = list(map(get_contiguous_orthologues, genome1_ortho_combinations))

    orthologue_number = len(orthologues)
    contiguous_orthologues_number = len(contiguous_orthologues)+1
    if orthologue_number:
        gene_order_conservation_score = float(contiguous_orthologues_number) / orthologue_number
        return gene_order_conservation_score, contiguous_orthologues_number, orthologue_number, orthologues
    return 0.0, contiguous_orthologues_number, orthologue_number, orthologues

def update_average_normalized_bit_score(fig_object_dict):

    for fig in fig_object_dict:
        for query in fig_object_dict[fig].peg_candidates:
            if query in fig_object_dict:
                average_normalized_bit_score = get_average_normalized_bit_score(fig_object_dict[fig],
                                                                                fig_object_dict[query])

                # add average_normalized_bit_score in fig_object.peg_candidates
                fig_object_dict[fig].peg_candidates[query]['average_normalized_bit_score'] = average_normalized_bit_score
    return fig_object_dict

def get_fig_BBHs(fig_object1, fig_object_dict):

    for query in fig_object1.peg_candidates:
        if query in fig_object_dict and query != fig_object1.ID:
            bol = are_bidirectional_best_hit(fig_object1,
                                             fig_object_dict[query],
                                             fig_object_dict)
            fig_object1.peg_candidates[query]['bidirectional_best_hit'] = bol
            if bol:
                yield query

def get_fig_triangle(fig_object1, ref_fig_object_dict, fig_object_dict):

    genome1 = get_orgid_from_fig(fig_object1.ID)

    query_list = list([query for query in fig_object1.peg_candidates.keys()
                       if get_orgid_from_fig(query) != genome1])
    for i in range(len(query_list)):
        for j in range(len(query_list)):
            if i<j:
                bol = are_TRIANGLE(fig_object1,
                                   ref_fig_object_dict[query_list[i]],
                                   ref_fig_object_dict[query_list[j]],
                                   fig_object_dict)
                fig_object1.peg_candidates[query_list[i]]['triangle'] = bol
                fig_object1.peg_candidates[query_list[j]]['triangle'] = bol
                if bol:
                    yield (query_list[i], query_list[j])

def get_fig_orthologues(fig_object, fig_object_dict, orthologuous_thresholds):

    genome1 = get_orgid_from_fig(fig_object.ID)
    for query in fig_object.peg_candidates:
        if query in fig_object_dict and get_orgid_from_fig(query) != genome1:
            orthologuous, aa_identity, length_coverage = is_orthologuous(fig_object,
                                                                         fig_object_dict[query],
                                                                         orthologuous_thresholds)
            if orthologuous:
                yield query

def get_gsp_confidence_score(fig_object_dict):

    gsps = [fig_object.ID for fig_object in fig_object_dict.values() if fig_object.ev]
    gsp_confidence_score = len(gsps) / len(fig_object_dict)
    return gsps, gsp_confidence_score

def annotation_confidence_score(fig_object1, fig_object_dict, ref_fig_object_dict, orthologuous_thresholds):

    orthologues = [elt for elt in get_fig_orthologues(fig_object1, fig_object_dict, orthologuous_thresholds)]
    BBHs = [elt for elt in get_fig_BBHs(fig_object1, fig_object_dict)]
    triangle = [elt for elt in get_fig_triangle(fig_object1, ref_fig_object_dict, fig_object_dict)]

    print(fig_object1.ID)
    ortho_BBHs = set(BBHs).intersection(orthologues)
    print(ortho_BBHs)
    pprint([ref_fig_object_dict[fig].peg for fig in ortho_BBHs if fig in ref_fig_object_dict])


def write_fasta_sequences(list_figs, f_out):
    assert type(list_figs) == list

    with open(f_out, 'w') as f1:
        for fig in list_figs:
            fig_info = query_fasta_sequences(fig)
            fasta_info = '> {0} {1} []\n{2}\n'.format(fig_info['Fig'], fig_info['Prot_name'], fig_info['Seq'])
            f1.write(fasta_info)
    return f_out


if __name__ == "__main__":

    candidate_fig_object_dict = pickle.load(open('Inputs/pickle/candidate_fig_object_dict.pkl', "rb"))
    ref_fig_object_dict = pickle.load(open('Inputs/pickle/ref_fig_network1_object_dict.pkl', "rb"))



