# SemiAutomaticPubSeedAnnotation

## Aims
1. Automatically annotate newly added genomes based on annotation data produced by the
laboratory
2. Set-up a high-quality automatic/semi-automatic annotation process
3. Rely on subsystems approaches to annotate sets of logically related functional roles
4. Facilitate comparative genomic analysis of annotations produced by the laboratory

## Strategies
1. Use subsystems to make easier the annotation process
2. Prioritize the Gold Standard Proteins for predictions
3. Construct an incremental process
4. Prioritize clustering-based annotation transfer over homology-based annotation
transfer
5. Respect the annotations standards of PEGs from GeneOntology + UniprotKB
6. Build small-scale annotation process from reduced reference dataset

## Requirements
### Tools + DBs
- Usearch
- Access to Gold Standard Proteins containing databases : UniProtKB, EcoCyc, MetaCyc, BRENDA, MicrobesOnline
- Access to PubSeed database
### Raw data
- Annotation spreadsheets from PubSeed subsystem of interest (.tsv) (i.e. http://pubseed.theseed.org//SubsysEditor.cgi?page=ShowSpreadsheet&subsystem=Propionate_production_HGM)
- Table of enzymes necessary for the functioning of a metabolic pathway (.csv)
- Table of enzyme subunits necessary for the functioning of an enzyme (.csv)
- List of Protein-encoding-genes to annotate (refered as PEGs of interest)

## SemiAutomaticPubSeedAnnotation Pipeline
### Overview
1. Find missing and mis-annotated Protein-Encoding Genes in the annotation spreadsheets from PubSeed subsystem
2. From the given list of Protein-Encoding-Genes (PEGs) to annotate, collect, from Gold Standard Proteins containing databases, Gold Standard Proteins refering to PEGs and homologs
3. From the given list of Protein-Encoding-Genes (PEGs) to annotate, collect PEGs from Annotation spreadsheets from PubSeed subsystem
4. Merge collected data from (2., 3.), refer them as reference PEGs, build collection of reference clusters from these PEGs
5. Check for the consistency of the annotations regarding the annotations standards
6. With given genome to annotate, identify candidate clusters that might correspond to PEGs of interest
7. Evaluate candidate clusters for their similarity with reference clusters built in (4.)
8. If the similarity of a candidate cluster with a refernece cluster is relevant, annotate all PEGs of the candidate cluster with the reference cluster PEG annotations
9. Hightlight pegs for which prediction could not distinguish between different reference proteins for manual annotation

## PubSeed automatic curation of Gold Standard Proteins
### Aim
(1) From set of genes, automatically collect their matching Gold Standard Proteins and all their homologs from Gold Standard Proteins containing databases (UniProtKB, EcoCyc, MetaCyc, BRENDA, MicrobesOnline) \
(2) Curate PubSeed database for pubmed_IDs that experimentally characterized the function of the PEG and check the consistency of the PEG function regarding annotation standards (Gene Ontology + UniprotKB) \
### Usage
PaperBlast.py [-h] -l LIST -d DEPTH [-p PATHWAY] \
LIST, help='str(list of proteins of interest)', type=str, required=True) i.e. "EG10585, Moth_1192, Moth_1193" \
DEPTH, help='depth defines the number of times newly found Gold Standard Proteins are re-used to find new homologs', type=int, required=True) \
PATHWAY, help='name of pathway', type=str) i.e. "propionate production" \


