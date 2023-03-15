#!/usr/bin/env python3

import argparse
from collections import Counter
from glob import glob
import json
import logging
import multiprocessing
import pathlib
import subprocess

from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np

from utils import fasta_iter

logger = multiprocessing.log_to_stderr()
logger.setLevel(logging.INFO)

class Protein():

    STANDARD_AMINO_ACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

    NH2O_RQEC = {
        'A' : 0.369, 'C' : -0.025, 'D' : -0.122, 
        'E' : -0.107, 'F' : -2.568, 'G' : 0.478, 
        'H' : -1.825, 'I' : 0.660, 'K' : 0.763, 
        'L' : 0.660, 'M' : 0.046, 'N' : -0.122, 
        'P' : -0.354, 'Q' : -0.107, 'R' : 0.072, 
        'S' : 0.575, 'T' : 0.569, 'V' : 0.522, 
        'W' : -4.087, 'Y' : -2.499
        }

    WEIGHTED_ZC = {
        'A': 0, 'C': 2.0, 'D': 4,
        'E': 2.0, 'F': -4.0, 'G': 2,
        'H': 4.0, 'I': -6, 'K': 4.0,
        'L': -6, 'M': -1.6, 'N': 4,
        'P': -2.0, 'Q': 2.0, 'R': 2.0,
        'S': 1.98, 'T': 0, 'V': -4.0,
        'W': -2.0, 'Y': -2.0
        }

    def __init__(self, protein_sequence : str):
        """
        :param protein_sequence: str
            Amino acid sequence of one protein
        """
        self.sequence =  self._format_protein_sequence(protein_sequence)
        self.length = len(self.sequence)

    def _format_protein_sequence(self, protein_sequence : str) -> str:
        """Returns a formatted amino acid sequence"""
        return ''.join([aa for aa in protein_sequence.strip().upper() if aa in self.STANDARD_AMINO_ACIDS])

    def _sequence_weighted_average(self, dict_ : dict):
        """Sums values in a dictionary and divides by sequence length"""
        return sum([dict_[s] for s in self.sequence]) / self.length

    def fraction_thermostable_ivywrel(self) -> float:
        """
        Thermostable residues reported by: 
        # https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030005
        """
        if self.length > 0:
            thermostable_length = len([aa for aa in self.sequence if aa in {'I', 'V', 'Y', 'W', 'R', 'E', 'L'}])
            return thermostable_length / self.length
        else:
            return np.nan
        
    def isoelectric_point(self) -> float:
        """Compute the isoelectric point (pI) of the protein"""
        if self.length > 0:
            # to-do: remove unnecessary Biopython dependency
            return IsoelectricPoint(self.sequence).pi()
        else:
            return np.nan
        
    def gravy(self):
        """Compute the Grand Average of Hydropathy (GRAVY)"""
        if self.length > 0:
            # to-do: remove unnecessary Biopython dependency
            return ProteinAnalysis(self.sequence).gravy()
        else:
            return np.nan
        
    def zc(self) -> float:
        """
        Computes average carbon oxidation state (Zc) of a 
        protein based on a dictionary of amino acids.
        """
        return self._sequence_weighted_average(dict_=self.WEIGHTED_ZC)

    def nh2o(self) -> float:
        """
        Computes  stoichiometric hydration state (nH2O) of a 
        protein based on a dictionary of amino acids.
        """
        return self._sequence_weighted_average(dict_=self.NH2O_RQEC)

class Genome():

    def __init__(self, protein_fasta_filepath : pathlib.Path):
        self.faa_filepath = str(protein_fasta_filepath)
        self.prefix = '.'.join(self.faa_filepath.split('/')[-1].split('.')[:-1])


    def _nonnull(list_ : list):
        """Returns array without NaN values"""
        X = np.array(list_)
        return X[~X.isnull()]

    def _length_weighted(statistics_lists : dict, property : str) -> float:
        sum_ = np.mean(_nonnull([f * l for f, l in zip(statistics_lists[property], statistics_lists['length'])]))
        sum_weights = np.sum(_nonnull(statistics_lists['length']))
        if sum_weights > 0:
            return sum_ / sum_weights
        else:
            return np.nan

    def _bin_midpoints(data : np.array, bins: np.array):
        """Returns counts per bin keyed by bin midpoint"""
        bin_counts = []
        bin_midpoints = []
        counts = Counter(np.digitize(data, bins))
        for bin_idx in range(len(bins) - 1):
            bin_counts.append(counts.get(bin_idx, 0))
            bin_midpoints.append(np.mean([bins[bin_idx], bins[bin_idx + 1]]))
        return dict(zip(bin_midpoints, bin_counts))

    def measure_protein_properties(self):
        """
        Returns a dictionary of properties for each protein.

        Association of properties to every protein allows statistics 
        to be computed jointly with multiple values, such as weighting
        a statistic by protein length.
        
        """
        protein_statistics_dict = {}
        with open(self.faa_filepath, 'r') as fh:
            for header, sequence in fasta_iter(fh):
                protein_id = header.split(' ')[0]
                protein_calc = Protein(sequence)
                # to-do: add other statistics below
                protein_statistics_dict[protein_id] = {
                    'length' : protein_calc.length,
                    'pi' : protein_calc.isoelectric_point(),
                    'gravy' : protein_calc.gravy(),
                    'zc' : protein_calc.zc(),
                    'nh2o' : protein_calc.nh2o(),
                    'f_ivywrel' : protein_calc.fraction_thermostable_ivywrel(),
                }
                
        return protein_statistics_dict

    def collect_genomic_statistics(self) -> dict:
        """
        Returns a dictionary of genome-wide statistics, based on 
        measurements, to be used for downstream analyses
        """
        logger.info("{}: Collecting protein statistics".format(self.prefix))
        try:
            genomic_properties = {}         
            protein_statistics_dict = self.measure_protein_properties()
            # may later measure other stats

            # reformat
            statistics_lists = defaultdict(list)
            for protein, stats in protein_statistics_dict.items():
                for key, value in stats.items():
                    statistics_lists[key].append(value)

            # means
            genomic_properties['mean_protein_length'] = np.mean(_nonnull(statistics_lists['length']))
            genomic_properties['mean_pi'] = np.mean(pis)
            genomic_properties['mean_gravy'] = np.mean(_nonnull(statistics_lists['gravy']))
            genomic_properties['mean_zc'] = np.mean(_nonnull(statistics_lists['zc']))
            genomic_properties['mean_nh2o'] = np.mean(_nonnull(statistics_lists['nh2o']))
            genomic_properties['mean_f_ivywrel'] = np.mean(_nonnull(statistics_lists['f_ivywrel']))

            # weighted means
            genomic_properties['weighted_mean_f_ivywrel'] = _length_weighted(statistics_lists, 'f_ivywrel')
            genomic_properties['weighted_mean_zc'] = _length_weighted(statistics_lists, 'zc')
            genomic_properties['weighted_mean_nh2o'] = _length_weighted(statistics_lists, 'nh2o')
            genomic_properties['weighted_mean_gravy'] = _length_weighted(statistics_lists, 'gravy')

            # distributions
            pis = _nonnull(statistics_lists['pi'])
            genomic_properties['histogram_pi'] = _bin_midpoints(pi, bins=np.linspace(0, 14, 141))
            genomic_properties['ratio_acidic_pis'] = len(pis[pis < 7]) / len(pis[pis >= 7])
            
            return genomic_properties
        except:
            return {}

def _mapping_wrapper(path):
    prefix = str(path).split('/')[-1].replace('_protein.faa', '')
    return prefix, Genome(path).collect_genomic_statistics()

def _format_pathlist_input(txt_file : str) -> list:
    pathlist = []
    with open(txt_file) as fh:
        for line in fh.readlines():
            path = line.strip()
            if pathlib.Path(path).is_file():
                pathlist.append(path)
    return list(set(pathlist))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                    prog='GenomicProperties',
                    description='Computes statistics from genomes supplied in FASTA format'
                    )
    
    parser.add_argument('-faas', help='Path to subdirectory with proteins as amino acids in FASTA format')
    parser.add_argument('-p', default=4, help='Number of parallel processes', required=False)
    parser.add_argument('-o', '--output', default='genomic_properties.json', help='Output file name, default fasta_prefix.json')

    args = parser.parse_args()
    
    pathlist = [path for path in pathlib.Path(args.faas).rglob('*.faa')] # _format_pathlist_input(str(args.genomes))
    output_json=str(args.output)
    workers = int(args.p)

    if workers is None:
        workers = multiprocessing.cpu_count() - 1

    logger.info("Measuring {} genomes with {} CPUs".format(len(pathlist), workers))
    filepath_gen = (pathlib.Path(path) for path in pathlist)
    with multiprocessing.Pool(workers) as p:
        pipeline_gen = p.map(_mapping_wrapper, filepath_gen)
        genomic_properties = dict(pipeline_gen)
    logger.info("Measured {} genomes".format(len(genomic_properties.keys())))
    
    json.dump(genomic_properties, open(output_json, 'w'))