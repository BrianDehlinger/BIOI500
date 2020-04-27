import json
import os
import pandas as pd
from collections import namedtuple
import re
import copy
from functools import reduce
import argparse


def add_gene(accession, gene_name, resistance_type, tool, gene_holder):
    if 'quinolone' in resistance_type.lower() or 'ciprofloxacin' in resistance_type.lower():
        tool[accession]['quinolone'].add(gene_name)
        gene_holder['quinolone'].add(gene_name)
    elif 'sulfonamide' in resistance_type.lower() or 'sulphonamide' in resistance_type.lower():
        tool[accession]['sulfonamide'].add(gene_name)
        gene_holder['sulfonamide'].add(gene_name)
    elif 'trimethoprim' in resistance_type.lower():
        tool[accession]['trimethoprim'].add(gene_name)
        gene_holder['trimethoprim'].add(gene_name)
    elif 'fosfomycin' in resistance_type.lower():
        tool[accession]['fosfomycin'].add(gene_name)
        gene_holder['fosfomycin'].add(gene_name)
    else:
        tool[accession]['other'].add(gene_name)
        gene_holder['other'].add(gene_name)

Resistance = namedtuple('Resistance', ('resistance_gene', 'predicted_phenotype', 'type'))

gene_map = {'quinolone': set(), 'sulfonamide': set(), 'trimethoprim': set(), 'fosfomycin': set(), 'other': set()}
amr_finder_genes = {}
amr = copy.deepcopy(gene_map)
ariba = copy.deepcopy(gene_map)
resfinder = copy.deepcopy(gene_map)
ariba_genes = {}
resfinder_genes = {}

parser = argparse.ArgumentParser()
parser.add_argument('--genomes', help="directory where genomes are located", required = True)
parser.add_argument('--prefix', help='prefix for resfinders data', required=False)
parser.add_argument('-noi', help='consider I to be susceptbile', required=False, action='store_true')
parser.add_argument('-allsnps', help='consider even unknown snps for ARIBA', required=False, action='store_true')
args = parser.parse_args()
genome_file_paths = os.listdir(args.genomes)

phenotypes = {}
actual_ecoli_phenotype_df = pd.read_csv("merged_ecoli_data")
accession_map = {}
reverse_map = {}
for index, row in actual_ecoli_phenotype_df.iterrows():
    accession_map[str(row['Strain'])] = row['RefSeq Assembly Accession #']
    reverse_map[str(row['RefSeq Assembly Accession #'])] = str(row['Strain'])
    amr_finder_genes[row['RefSeq Assembly Accession #']] = copy.deepcopy(gene_map)
    ariba_genes[row['RefSeq Assembly Accession #']] = copy.deepcopy(gene_map)
    resfinder_genes[row['RefSeq Assembly Accession #']] = copy.deepcopy(gene_map)

actual_ecoli_phenotype_df = actual_ecoli_phenotype_df.drop(['AMC', 'CPD', 'SRA Accession #', 'Strain', 'Patient Symptom'], axis=1)
actual_ecoli_phenotype_df.set_index('RefSeq Assembly Accession #', inplace=True)
ariba_ecoli_data = pd.read_csv('e_coli_resfinder.csv')
ariba_parameter = pd.read_csv('e_coli_parameter.csv')
ariba_ecoli_data = ariba_ecoli_data[['name', 'sul1.match', 'sul2.match', 'qnrB-.match', 'dfrA8.match', 'dfrA2-.match', 'dfrA14.match', 'dfrA-.match', 'aac_6___Ib_cr+.match']]
ariba_parameter = ariba_parameter[['name', 'sul1.match', 'sul2.match', 'dfrA8.match', 'dfrA2-.match', 'dfrA14.match', 'dfrA-.match']]
ariba_ecoli_data['name'] = ariba_ecoli_data['name'].apply(lambda string: string.split('/')[0])
ariba_parameter['name'] = ariba_parameter['name'].apply(lambda string: string.split('/')[0])
if args.allsnps:
    ariba_ecoli_snps = pd.read_csv('e_coli_all_snps.csv')
else:
    ariba_ecoli_snps = pd.read_csv('e_coli_snps.csv')
columns = list(ariba_ecoli_snps.columns)
ariba_ecoli_snps['name'] = ariba_ecoli_snps['name'].apply(lambda string: string.split('/')[0])
ariba_ecoli_data = pd.merge(ariba_ecoli_snps, ariba_ecoli_data, left_on='name', right_on='name', how='inner')
ariba_parameter = pd.merge(ariba_parameter, ariba_ecoli_snps, left_on='name', right_on='name', how='inner')
ariba_ecoli_data['accession'] = ariba_ecoli_data['name'].apply(lambda string: accession_map[string])
ariba_parameter['accession'] = ariba_parameter['name'].apply(lambda string: accession_map[string])

for index, row in ariba_ecoli_data.iterrows():
    trimethoprim = [row['dfrA2-.match'], row['dfrA14.match'], row['dfrA2-.match'], row['dfrA-.match']]
    sulfonamide = [row['sul2.match'], row['sul1.match']]
    quinolone = [row['qnrB-.match'], row['aac_6___Ib_cr+.match']] + [row[item] for item in columns]

for index, row in ariba_parameter.iterrows():
    trimethoprim = [row['dfrA2-.match'], row['dfrA14.match'], row['dfrA2-.match'], row['dfrA-.match']]
    sulfonamide = [row['sul2.match'], row['sul1.match']]
    quinolone = [row[item] for item in columns]

if args.noi:

    def add_actual_phenotype_ecoli(accession, actual_dictionary):
        resistance_set = actual_dictionary
        row = actual_ecoli_phenotype_df.loc[accession]
        fosfomycin_resistance = row.FOF
        quinolone_resistance = row.CIP
        sulfonamide_and_trimethoprim_resistance = row.SXT
        if sulfonamide_and_trimethoprim_resistance != 'S' and sulfonamide_and_trimethoprim_resistance != 'I':
            resistance_set.add('sulfonamide resistance')
            resistance_set.add('trimethoprim resistance')
        if quinolone_resistance != 'S' and quinolone_resistance != 'I':
            resistance_set.add('quinolone resistance')
        if fosfomycin_resistance != 'S' and fosfomycin_resistance != 'I':
            resistance_set.add('fosfomycin resistance')
else:

    def add_actual_phenotype_ecoli(accession, actual_dictionary):
        resistance_set = actual_dictionary
        row = actual_ecoli_phenotype_df.loc[accession]
        fosfomycin_resistance = row.FOF
        quinolone_resistance = row.CIP
        sulfonamide_and_trimethoprim_resistance = row.SXT
        if sulfonamide_and_trimethoprim_resistance != 'S':
            resistance_set.add('sulfonamide resistance')
            resistance_set.add('trimethoprim resistance')
        if quinolone_resistance != 'S':
            resistance_set.add('quinolone resistance')
        if fosfomycin_resistance != 'S':
            resistance_set.add('fosfomycin resistance')

for genome in genome_file_paths:
    amrfinder_resistance_phenotypes = set() 
    resfinder_phenotypes = set()
    phenotypes[genome] = {'amrfinder': set(), 'resfinder': set(), 'actual': set()}
    add_actual_phenotype_ecoli(genome, phenotypes[genome]['actual'])
    # Resfinder -> Finding ARGs

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
                   add_gene(genome, resistance_payload['resistance_gene'], predicted_phenotype, resfinder_genes, resfinder)

    # Pointfinder -> Point mutations that confer resistance
    if args.prefix is None:
        df = pd.read_csv(f'{args.genomes}/{genome}/resfinder/GCF_blastn_results.tsv', '\t')
    if args.prefix is not None:
        df = pd.read_csv(f'{args.genomes}/{genome}/resfinder/{genome}_blastn_results.tsv', '\t')
    map = {'Nalidixic acid': 'Quinolone resistance', 'Ciprofloxacin': 'Quinolone resistance'}
    if df.empty:
        pass
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
                    print(resistance)
                    list_of_resfinder_resistances.append(resistance)
                    resfinder_phenotypes.add(phenotype)
                    add_gene(genome, gene, phenotype, resfinder_genes, resfinder)
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
                add_gene(genome, gene, phenotype, resfinder_genes, resfinder)
 
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
            print(f'RESISTANCE: {resistance} {reverse_map[genome]}')
            list_of_amr_resistances.append(resistance)
            amrfinder_resistance_phenotypes.add(antibiotic_class.lower() + ' resistance')
            add_gene(genome, gene, antibiotic_class.lower(), amr_finder_genes, amr)

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

def _analyze_trimethoprim_and_sulfonamide(actual, predicted):
    if 'trimethoprim resistance' in actual and 'sulfonamide resistance' in actual:
        if 'trimethoprim resistance' in predicted and 'sulfonamide resistance' in predicted:
            return 'correct_in'
        else:
            return 'false_negative'
    else:
        if 'trimethoprim_resistance' not in predicted and 'sulfonamide resistance' not in predicted:
            return 'correct_not_in'
        else:
            return 'false_positive'

def _analyze_quinolone(actual, predicted):
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

def _analyze_fosfomycin(actual, predicted):
    if 'fosfomycin resistance' in actual:
        if 'fosfomycin resistance' in predicted:
            return 'correct_in'
        else:
            return 'false_negative'
    else:
        if 'fosfomycin resistance' not in predicted:
            return 'correct_not_in'
        else:
            return 'false_positive'

def analyze_amr_e_coli(actual, predicted):
    fosfomycin_bundle = _analyze_fosfomycin(actual, predicted)
    quinolone_bundle = _analyze_quinolone(actual, predicted)
    trimethoprim_and_sulfonamide_bundle = _analyze_trimethoprim_and_sulfonamide(actual, predicted)
    return fosfomycin_bundle, quinolone_bundle, trimethoprim_and_sulfonamide_bundle

base_drug_map = {'false_negative': 0, 'false_positive': 0, 'correct_in': 0, 'correct_not_in': 0}
drug_map = {'trimethoprim_and_sulfonamide': copy.deepcopy(base_drug_map), 'quinolone': copy.deepcopy(base_drug_map), 'fosfomycin': copy.deepcopy(base_drug_map)}
stats = {'resfinders': copy.deepcopy(drug_map), 'amrfinder': copy.deepcopy(drug_map), 'ariba': copy.deepcopy(drug_map), 'ariba_parameter': copy.deepcopy(drug_map)}
for accession, phenotype_data in phenotypes.items():
    actual = phenotype_data['actual']
    fosfomycin_res, quinolone_res, trimethoprim_and_sulfonamide_res= analyze_amr_e_coli(actual, phenotype_data['resfinder'])
    print(reverse_map[accession], phenotype_data['amrfinder'])
    fosfomycin, quinolone, trimethoprim_and_sulfonamide = analyze_amr_e_coli(actual, phenotype_data['amrfinder'])
    stats['resfinders']['fosfomycin'][fosfomycin_res] += 1
    stats['resfinders']['quinolone'][quinolone_res] += 1
    stats['resfinders']['trimethoprim_and_sulfonamide'][trimethoprim_and_sulfonamide_res] += 1
    stats['amrfinder']['fosfomycin'][fosfomycin] += 1
    stats['amrfinder']['quinolone'][quinolone] += 1
    stats['amrfinder']['trimethoprim_and_sulfonamide'][trimethoprim_and_sulfonamide] += 1

ariba_ecoli_data = ariba_ecoli_data.loc[ariba_ecoli_data['accession'].isin(genome_file_paths)]
ariba_parameter = ariba_parameter.loc[ariba_parameter['accession'].isin(genome_file_paths)]
for index, row in ariba_ecoli_data.iterrows():
    trimethoprim = any([item == 'yes' for item in [row['dfrA2-.match'], row['dfrA14.match'], row['dfrA2-.match'], row['dfrA-.match']]])
    sulfonamide = any([item == 'yes' for item in [row['sul2.match'], row['sul1.match']]])
    quinolone = [row['qnrB-.match'], row['aac_6___Ib_cr+.match']] + [row[item] for item in columns]
    quinolone = any([item == 'yes' for item in quinolone])
    resistances = set()
    if quinolone:
        resistances.add('quinolone resistance')
    if sulfonamide:
        resistances.add('sulfonamide resistance')
    if trimethoprim:
        resistances.add('trimethoprim resistance')
    fosfomycin, quinolone, trimethoprim_and_sulfonamide = analyze_amr_e_coli(phenotypes[row.accession]['actual'], resistances)
    stats['ariba']['fosfomycin'][fosfomycin] += 1
    stats['ariba']['quinolone'][quinolone] += 1
    stats['ariba']['trimethoprim_and_sulfonamide'][trimethoprim_and_sulfonamide] += 1


for index, row in ariba_parameter.iterrows():
    print(row.accession)
    trimethoprim = any([item == 'yes' for item in [row['dfrA2-.match'], row['dfrA14.match'], row['dfrA2-.match'], row['dfrA-.match']]])
    sulfonamide = any([item == 'yes' for item in [row['sul2.match'], row['sul1.match']]])
    quinolone = [row[item] for item in columns]
    quinolone = any([item == 'yes' for item in quinolone])
    resistances = set()
    if quinolone:
        resistances.add('quinolone resistance')
    if sulfonamide:
        resistances.add('sulfonamide resistance')
    if trimethoprim:
        resistances.add('trimethoprim resistance')
    fosfomycin, quinolone, trimethoprim_and_sulfonamide = analyze_amr_e_coli(phenotypes[row.accession]['actual'], resistances)
    stats['ariba_parameter']['fosfomycin'][fosfomycin] += 1
    stats['ariba_parameter']['quinolone'][quinolone] += 1
    stats['ariba_parameter']['trimethoprim_and_sulfonamide'][trimethoprim_and_sulfonamide] += 1

print(stats)

# Checking for identical called genes for each class of genes
for accession in amr_finder_genes:
    quinolones = sorted(list(amr_finder_genes[accession]['quinolone']))
    print(f"Record {accession}:\n\n")
    print(f"There are {len(quinolones)} quinolone genes for AMRFinderPlus!")
    resfinder_quinolones = sorted([item.replace(' p.', '_') for item in list(resfinder_genes[accession]['quinolone'])])
    resfinder_quinolones = sorted([item.replace('aac(6\')-Ib-cr', 'aac(6\')-Ib-cr5') for item in resfinder_quinolones])
    if (resfinder_quinolones != quinolones):
        print(resfinder_quinolones, quinolones)
    else:
        print("True")
    trimethoprim = sorted(list(amr_finder_genes[accession]['trimethoprim']))
    resfinder_trimethoprim = sorted(resfinder_genes[accession]['trimethoprim'])
    sulfonamide = sorted(list(amr_finder_genes[accession]['sulfonamide']))
    resfinder_sulfonamide = sorted(resfinder_genes[accession]['sulfonamide'])
    if resfinder_trimethoprim != trimethoprim:
        print(resfinder_trimethoprim, trimethoprim)
    if resfinder_sulfonamide != sulfonamide:
        print(resfinder_sulfonamide, sulfonamide)
    fosfomycin = sorted(list(amr_finder_genes[accession]['fosfomycin']))
    resfinder_fosfomycin = sorted(resfinder_genes[accession]['fosfomycin'])
    if (fosfomycin != resfinder_fosfomycin):
        print("False")
        print(resfinder_fosfomycin, fosfomycin)
    else:
        print("True")


print(amr)
print(resfinder)

