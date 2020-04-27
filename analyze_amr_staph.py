import json
import os
import pandas as pd
from collections import namedtuple
import re
import copy
from functools import reduce
import argparse

Resistance = namedtuple('Resistance', ('resistance_gene', 'predicted_phenotype', 'type'))

parser = argparse.ArgumentParser()
parser.add_argument('--genomes', help='directory where genomes are located', required=True)
parser.add_argument('-allsnps', help='consider even unknown snps for ARIBA', required=False, action='store_true')
parser.add_argument('--mic', help='mic to consider resistant', required=True)
parser.add_argument('-ariba_only_two_snps', help='fair comparison by using only 2SNPS present in other databases for ARIBA', action='store_true')
args = parser.parse_args()
genome_file_paths = os.listdir(args.genomes)
genome_file_paths.remove('AA104')
genome_file_paths.remove('AA95')

def extract_AA(string):
    item = string.split('/')
    if 'ariba_snp' in item[1]:
        return string.split('/')[0]
    elif 'snp' in item[0]:
        return string.split('_')[0]

phenotypes = {}
actual_staph_phenotypes = pd.read_csv("staph_phenotype.csv")
actual_staph_phenotypes = actual_staph_phenotypes[actual_staph_phenotypes['accession'] != 'AA104']
actual_staph_phenotypes = actual_staph_phenotypes[actual_staph_phenotypes['accession'] != 'AA95']
actual_staph_phenotypes.set_index('accession', inplace=True)
ariba_staph_data = pd.read_csv("staph_snps.csv")
if args.allsnps:
    ariba_staph_data = pd.read_csv("staph_all_snps")
if args.ariba_only_two_snps:
    ariba_staph_data = ariba_staph_data[['name', 'gyrA.S84L', 'grlA.S80F']]
ariba_staph_data['old_name'] = ariba_staph_data['name']
ariba_staph_data['name'] = ariba_staph_data['name'].apply(lambda string: extract_AA(string))
ariba_staph_data['accession'] = ariba_staph_data['name']
ariba_staph_data = ariba_staph_data[ariba_staph_data['accession'] != 'AA104']
ariba_staph_data = ariba_staph_data[ariba_staph_data['accession'] != 'AA95']
columns = list(ariba_staph_data.columns)
for index, row in ariba_staph_data.iterrows():
    quinolone = [row[item] for item in columns]

def add_actual_phenotype(accession, actual_dictionary):
    resistance_set = actual_dictionary
    row = actual_staph_phenotypes.loc[accession]
    quinolone_resistance = row.CIP
    print(row)
    print(quinolone_resistance)
    if int(quinolone_resistance) >= int(args.mic):
        resistance_set.add('quinolone resistance')
    else:
        print("NOT")

for genome in genome_file_paths:
    amrfinder_resistance_phenotypes = set()
    resfinder_phenotypes = set()
    phenotypes[genome] = {'amrfinder': set(), 'resfinder': set(), 'actual': set()}
    add_actual_phenotype(genome, phenotypes[genome]['actual'])

    list_of_resfinder_resistances = []
    list_of_amr_resistances = []

    with open(f'{args.genomes}/{genome}/resfinder/data_resfinder.json', 'r') as resfinder_file:
        data = resfinder_file.read()
        resfinder_json = json.loads(data)['resfinder']

        results_json = resfinder_json['results']
        for antibiotic_type, result_json in results_json.items():
            if result_json[antibiotic_type.lower()] == 'No hit found':
                pass
            else:
                for hit, resistance_payload in result_json[antibiotic_type.lower()].items():
                   predicted_phenotype = resistance_payload['predicted_phenotype']
                   if 'Warning' in predicted_phenotype:
                       predicted_phenotype = antibiotic_type + ' resistance'
                   pattern = re.compile('.* [R, r]esistance')
                   predicted_phenotype = re.findall(pattern, predicted_phenotype)[0].lower()
                   resistance = Resistance(resistance_payload['resistance_gene'], predicted_phenotype, 'gene')
                   list_of_resfinder_resistances.append(resistance)
                   resfinder_phenotypes.add(predicted_phenotype)
    df = pd.read_csv(f'{args.genomes}/{genome}/resfinder/{genome}_blastn_results.tsv', '\t')
    map = {'Nalidixic acid': 'Quinolone resistance', 'Ciprofloxacin': 'Quinolone resistance'}
    if df.empty:
        print("DF is empty")
    else:
        for row in df.itertuples():
            gene = row.Mutation
            phenotype = row.Resistance
            if ',' in phenotype:
                phenotypes_local = phenotype.split(",")
                for phenotype in phenotypes_local:
                    if phenotype in map.keys():
                        phenotype = map[phenotype].lower()
                    pattern = re.compile('.* [R, r]esistance$')
                    phenotype = re.findall(pattern, phenotype)[0].lower()
                    resistance = Resistance(gene, phenotype, 'point')
                    list_of_resfinder_resistances.append(resistance)
                    resfinder_phenotypes.add(phenotype)
            else:
                resistance = Resistance(gene, phenotype.lower(), 'point')
                print(resistance)
                list_of_resfinder_resistances.append(resistance)
                pattern = re.compile('.* [R, r]esistance$')
                if 'resistance' not in phenotype.lower():
                    phenotype = phenotype + ' resistance'
                print(phenotype)
                phenotype = re.findall(pattern, phenotype)[0].lower()
                resfinder_phenotypes.add(phenotype)
    # AMRFinder -> Point mutations AND AMR genes from one command!
    try:
        df = pd.read_csv(f"{args.genomes}/{genome}/amrfinder/amrfinderoutput.txt", sep='\t')
    except:
        df = None
    if df is None:
        pass
    else:
        for row in df.itertuples():
            gene = row._6
            type = 'point' if row._10 == 'POINT' else 'gene'
            antibiotic_class = row.Class.lower()
            resistance = Resistance(gene, antibiotic_class.lower() + ' resistance', type = type)
            print(resistance)
            list_of_amr_resistances.append(resistance)
            amrfinder_resistance_phenotypes.add(antibiotic_class.lower() + ' resistance')

    if 'sulphonamide resistance' in amrfinder_resistance_phenotypes:
        amrfinder_resistance_phenotypes.remove('sulphonamide resistance')
        amrfinder_resistance_phenotypes.add('sulfonamide resistance')
    if 'sulphonamide resistance' in resfinder_phenotypes:
        resfinder_phenotypes.remove('sulphonamide resistance')
        resfinder_phenotypes.add('sulfonamide resistance')
    if 'quinolone/triclosan resistance' in amrfinder_resistance_phenotypes:
        amrfinder_resistance_phenotypes.remove('quinolone/triclosan resistance')
        amrfinder_resistance_phenotypes.add('quinolone resistance')
        amrfinder_resistance_phenotypes.add('triclosan resistance')
    if 'aminoglycoside/quinolone resistance' in amrfinder_resistance_phenotypes:
        amrfinder_resistance_phenotypes.remove('aminoglycoside/quinolone resistance')
        amrfinder_resistance_phenotypes.add('quinolone resistance')
        amrfinder_resistance_phenotypes.add('aminoglycoside resistance')
    if 'ciprofloxacin resistance' in amrfinder_resistance_phenotypes:
        amrfinder_resistance_phenotypes.add('quinolone resistance')
    if 'ciprofloxacin resistance' in resfinder_phenotypes:
        resfinder_phenotypes.add('quinolone resistance')
    phenotypes[genome]['resfinder'] = resfinder_phenotypes
    phenotypes[genome]['amrfinder'] = amrfinder_resistance_phenotypes

def _analyze_quinolone(actual, predicted):
    print(actual, predicted)
    if 'quinolone resistance' in actual:
        if 'quinolone resistance' in predicted:
            return 'correct_in'
        else:
            return 'false_negative'
    else:
        if 'quinolone resistance' not in predicted:
            return 'correct_not_in'
        else:
            return 'false_positive'

base_drug_map = {'false_negative': 0, 'false_positive': 0, 'correct_in': 0, 'correct_not_in': 0}
drug_map = {'quinolone': copy.deepcopy(base_drug_map)}
stats = {'resfinders': copy.deepcopy(drug_map), 'amrfinder': copy.deepcopy(drug_map), 'ariba': copy.deepcopy(drug_map)}

for accession, phenotype_data in phenotypes.items():
    actual = phenotype_data['actual']
    quinolone_resistance = _analyze_quinolone(actual, phenotype_data['resfinder'])
    amr_quinolone_resistance = _analyze_quinolone(actual, phenotype_data['amrfinder'])
    stats['resfinders']['quinolone'][quinolone_resistance] += 1
    stats['amrfinder']['quinolone'][quinolone_resistance] += 1



for index, row in ariba_staph_data.iterrows():
    quinolone = [row[item] for item in columns]
    quinolone = any([item == 'yes' for item in quinolone])
    resistances = set()
    if quinolone:
        resistances.add('quinolone resistance')
    quinolone = _analyze_quinolone(phenotypes[row.accession]['actual'], resistances)
    if quinolone == 'false_negative':
        print(f"FALSE NEGATIVE: {row.accession}")
    stats['ariba']['quinolone'][quinolone] += 1

print(stats)
