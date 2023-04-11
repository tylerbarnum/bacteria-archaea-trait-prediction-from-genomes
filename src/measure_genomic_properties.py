#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
from glob import glob
from itertools import chain, product
import json
import logging
import multiprocessing
from pathlib import Path
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
    
    HYDROPHOBICITY = {
        'A': 1.8, 'R': -4.5, 'N': -3.5,
        'D': -3.5, 'C': 2.5, 'Q': -3.5,
        'E': -3.5, 'G': -0.4, 'H': -3.2,
        'I': 4.5, 'L': 3.8, 'K': -3.9,
        'M': 1.9, 'F': 2.8, 'P': -1.6,
        'S': -0.8, 'T': -0.7, 'W': -0.9,
        'Y': -1.3,'V': 4.2
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

    def aa_frequencies(self) -> dict:
        """Returns count of every amino acid """
        return dict(Counter([aa for aa in self.sequence]))

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
            return np.mean([self.HYDROPHOBICITY[aa] for aa in self.sequence])
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

    def is_soluble_extracellular_heuristic(self) -> bool: 
        """
        Sec and Tat signal peptides are N-terminal sequences with a
        hydrophobic region. This function applies a heuristic for 
        the signal peptides: a region of high hydrophobicity at the 
        start protein sequence, optionally with an additional motif. 
        Sec signals are more hypophobic and Tat signals have a consensus
        motif of [S/T]RRxFLK.

        Compared to Signalp6, this heuristic is far, far faster. A
        downside is that, while it returns the majority of
        extracellular soluble proteins, it also returns about the
        same number of false positives. Therefore while this can be used to enrich
        for proteins that are extracellular and soluble, each protein only has
        a ~50% chance of being extracellular.
        """

        # Measure protein and N-terminal region
        max_signal_hydrophobicity = self.find_max_hydrophobicity(l_substring=8, l_peptide=35)
        max_tat_score = self.find_max_tat_motif_score(l_peptide=30)
        protein_hyprophob = self.gravy()
        
        # Score
        # Thresholds determined empirically but roughly
        threshold_hydrophob_soluble = 0 
        threshold_hydrophob_sec = 2.8
        threshold_hydrophob_tat = 1
        threshold_tat_score = 8.5
        if protein_hyprophob < threshold_hydrophob_soluble:
            if max_tat_score >= threshold_tat_score and max_signal_hydrophobicity >= threshold_hydrophob_tat:
                return True
            elif max_signal_hydrophobicity >= threshold_hydrophob_sec:
                return True
        return False

    def tat_motif_score(self, subsequence):
        """
        Heuristic to quickly score whether a TAT motif is present.
        Each position receives a score if an expected residue is
        present. Scores are based on consensus sequence motifs. Local
        alignment would be a more flexible improvement.
        """
        # Determined by eyeballing motifs
        position_scores = {
            0 : {'S' : 1.5, 'T' : 1.5},
            1 : {'R' : 4},
            2 : {'R' : 4},
            3 : {},
            4 : {'F' : 2, 'V' : 0.5, 'L' : 1},
            5 : {'L' : 3, 'V' : 0.5, 'I' : 0.5},
            6 : {'K' : 0.5}
        }
        score = 0
        for i, score_dict in position_scores.items():
            score += score_dict.get(subsequence[i], 0)
        return score

    def find_max_tat_motif_score(self, l_peptide : int = 30):
        """
        Scans the n-terminal region of a protein plausibly part of a 
        signal peptide and returns the maximum tat score within that region.
        """
        max_score = 0
        l_substring = 7
        for i in range(min([len(self.sequence) - l_substring + 1, l_peptide])):
            region = self.sequence[i:(i + l_substring)]
            region_score = self.tat_motif_score(region)
            if region_score >= max_score:
                max_score = region_score
        return max_score

    def find_max_hydrophobicity(self, l_substring : int = 10, l_peptide : int = 40):
        """
        Scans the n-terminal region of a protein plausibly part of a 
        signal peptide and returns the maximum hydrophobicity of a subsequence.
        """
        max_hydrophobicity = -999999 # nonsensical
        
        hydrophobicities = [self.HYDROPHOBICITY[aa] for aa in self.sequence]
        for i in range(min([len(hydrophobicities) - l_substring + 1, l_peptide])):
            region = hydrophobicities[i:(i + l_substring)]
            region_hydrophobicity = np.mean(region)
            if region_hydrophobicity >= max_hydrophobicity:
                max_hydrophobicity = region_hydrophobicity

        return max_hydrophobicity

class Nucleotide():

    def __init__(self, nucleotide : str):
        """
        :param nucleotide_sequence: str
            Nucleotide sequence
        :param k: int
            Length of k-mer to count
        """
        self.sequence =  nucleotide
        self.length = len(self.sequence)
    
    def count_canonical_kmers(self, k : int =2) -> dict:
        """
        Returns counts for canonical k-mers only, e.g. 'AA' records
        counts for both 'AA' and its reverse complement 'TT'.
        """
        kmers_count = self.count_kmers(k)
        
        canonical_kmers_dict = self.make_canonical_kmers_dict(k)
        canonical_kmers_count = defaultdict(int)
        for kmer, count in kmers_count.items():
            canonical_kmer = canonical_kmers_dict.get(kmer, None)
            canonical_kmers_count[canonical_kmer] += count
        return dict(canonical_kmers_count)

    def count_kmers(self, k : int =2) -> dict:
        """Returns count for each k-mer at specific k"""
        kmers_count = Counter([self.sequence[i:i+k] for i in range(len(self.sequence) - k + 1)])
        return dict(kmers_count)
         
    def make_canonical_kmers_dict(self, k : int=2):
        """
        Creates a dictionary that where a k-mer and its reverse
        complement k-mer are keys, and the k-mer is the value, for
        all k-mers at a given k. E.g.: {'AA' : 'TT', 'AA' : 'AA'}
        """
        canonical_kmers_dict = {}
        kmers_sorted = sorted([''.join(i) for i in product(['A', 'T', 'C', 'G'], repeat=k)])
        for kmer in kmers_sorted:
            if kmer not in canonical_kmers_dict.keys():
                canonical_kmers_dict[kmer] = kmer
                canonical_kmers_dict[self.reverse_complement(kmer)] = kmer
        return canonical_kmers_dict

    def reverse_complement(self, sequence):
        complement = {'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C', 'N' : 'N'}
        return ''.join([complement.get(nt, '') for nt in sequence[::-1]])



class Genome():

    def __init__(self, protein_faa_filepath : Path, fna_filepath : Path):
        self.faa_filepath = str(protein_faa_filepath)
        self.fna_filepath = str(fna_filepath)
        self.prefix = '.'.join(self.faa_filepath.split('/')[-1].split('.')[:-1])
        self.protein_statistics = None

    def measure_protein_statistics(self):
        """
        Returns a dictionary of properties for each protein.

        Association of properties to every protein allows statistics 
        to be computed jointly with multiple values, such as weighting
        a statistic by protein length.
        """
        if self.protein_statistics is None:
            self.protein_statistics = {}
            with open(self.faa_filepath, 'r') as fh:
                for header, sequence in fasta_iter(fh):
                    protein_id = header.split(' ')[0]
                    protein_calc = Protein(sequence)
                    # to-do: add other statistics below
                    self.protein_statistics[protein_id] = {
                        'length' : protein_calc.length,
                        'pi' : protein_calc.isoelectric_point(),
                        'gravy' : protein_calc.gravy(),
                        'zc' : protein_calc.zc(),
                        'nh2o' : protein_calc.nh2o(),
                        'aa_frequencies' : protein_calc.aa_frequencies(),
                        'f_ivywrel' : protein_calc.fraction_thermostable_ivywrel(),
                        'extra_soluble' : protein_calc.is_soluble_extracellular_heuristic(),
                    }
                    
        return self.protein_statistics

    def compute_proteome_statistics(self, subset_proteins : set=None) -> dict:
        """
        Returns a dictionary of genome-wide statistics, based on 
        measurements, to be used for downstream analyses
        """
        proteome_statistics = {}

        if subset_proteins:
            protein_statistics_dict = {k : self.measure_protein_statistics()[k] for k in subset_proteins}
        else:
            protein_statistics_dict = self.measure_protein_statistics()
        
        statistics_lists = defaultdict(list)
        for protein, stats in protein_statistics_dict.items():
            for key, value in stats.items():
                statistics_lists[key].append(value)

        # distributions / frequencies
        pis = self._nonnull(statistics_lists['pi'])
        proteome_statistics['histogram_pi'] = self._bin_midpoints(pis, bins=np.linspace(0, 14, 141))
        proteome_statistics['ratio_acidic_pis'] = len(pis[pis < 7]) / len(pis[pis >= 7])
        for aa, count in Counter(chain(*statistics_lists['aa_frequencies'])).items():
            proteome_statistics[f'aa_{aa.lower()}'] = count

        # means
        proteome_statistics['mean_protein_length'] = np.mean(self._nonnull(statistics_lists['length']))
        proteome_statistics['mean_pi'] = np.mean(pis)
        proteome_statistics['mean_gravy'] = np.mean(self._nonnull(statistics_lists['gravy']))
        proteome_statistics['mean_zc'] = np.mean(self._nonnull(statistics_lists['zc']))
        proteome_statistics['mean_nh2o'] = np.mean(self._nonnull(statistics_lists['nh2o']))
        proteome_statistics['mean_f_ivywrel'] = np.mean(self._nonnull(statistics_lists['f_ivywrel']))

        # weighted means
        proteome_statistics['weighted_mean_f_ivywrel'] = self._length_weighted(statistics_lists, 'f_ivywrel')
        proteome_statistics['weighted_mean_zc'] = self._length_weighted(statistics_lists, 'zc')
        proteome_statistics['weighted_mean_nh2o'] = self._length_weighted(statistics_lists, 'nh2o')
        proteome_statistics['weighted_mean_gravy'] = self._length_weighted(statistics_lists, 'gravy')

        return proteome_statistics

    def compute_genome_statistics(self):
        """
        Returns a dictionary of genome-wide statistics on nucleotide
        content, currently only k-mer counts
        """
        genome_statistics = {}
        with open(self.fna_filepath, 'r') as fh:
            genome = ''
            for header, sequence in fasta_iter(fh):
                genome += 'NN' + sequence
            nucleotide_calc = Nucleotide(genome)
            genome_statistics.update(nucleotide_calc.count_canonical_kmers(k=1))
            genome_statistics.update(nucleotide_calc.count_canonical_kmers(k=2))

        return genome_statistics

    def collect_genomic_statistics(self) -> dict:
        """
        Computes statistics about the proteome for each genome on:
        1. All proteins, keyed by 'all'
        2. Extracellular soluble proteins, keyed by 'extracellular_soluble'
        """
        
        self.genomic_statistics = {}

        logger.info("{}: Collecting protein statistics".format(self.prefix))
        self.genomic_statistics['all'] = self.compute_proteome_statistics()
        self.genomic_statistics['all'].update(self.compute_genome_statistics())

        # Extracellular and soluble proteins
        extra_sol_proteins = set()
        for protein, stats in self.measure_protein_statistics().items():
            if stats['extra_soluble'] == True:
                extra_sol_proteins.add(protein)
        self.genomic_statistics['extracellular_soluble'] = self.compute_proteome_statistics(subset_proteins=extra_sol_proteins)

        return self.genomic_statistics

    def _nonnull(self, list_ : list):
        """Returns array without NaN values"""
        X = np.array(list_)
        return X[~np.isnan(X)]

    def _length_weighted(self, statistics_lists : dict, property : str) -> float:
        sum_ = np.sum(self._nonnull([f * l for f, l in zip(statistics_lists[property], statistics_lists['length'])]))
        sum_weights = np.sum(self._nonnull(statistics_lists['length']))
        if sum_weights > 0:
            return sum_ / sum_weights
        else:
            return np.nan

    def _bin_midpoints(self, data : np.array, bins: np.array):
        """Returns counts per bin keyed by bin midpoint"""
        bin_counts = []
        bin_midpoints = []
        counts = Counter(np.digitize(data, bins))
        for bin_idx in range(len(bins) - 1):
            bin_counts.append(counts.get(bin_idx, 0))
            bin_midpoints.append(round(np.mean([bins[bin_idx], bins[bin_idx + 1]]), 3))
        return dict(zip(bin_midpoints, bin_counts))

def _mapping_wrapper(paths):
    faa_path, fna_path = paths
    prefix = str(faa_path).split('/')[-1].replace('_protein.faa', '')
    return prefix, Genome(protein_faa_filepath=faa_path, fna_filepath=fna_path).collect_genomic_statistics()

def _format_pathlist_input(txt_file : str) -> list:
    pathlist = []
    with open(txt_file) as fh:
        for line in fh.readlines():
            path = line.strip()
            if Path(path).is_file():
                pathlist.append(path)
    return list(set(pathlist))

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                    prog='GenomicProperties',
                    description='Computes statistics from genomes supplied in FASTA format'
                    )
    
    parser.add_argument('-faas', help='Path to subdirectory with proteins as amino acids in FASTA format, suffix .faa')
    parser.add_argument('-fnas', help='Path to subdirectory with proteins as nucleotides in FASTA format, suffix .fna')
    parser.add_argument('-p', default=4, help='Number of parallel processes', required=False)
    parser.add_argument('-o', '--output', default='genomic_properties.json', help='Output file name, default fasta_prefix.json')

    args = parser.parse_args()
    
    faa_pathlists = {str(path).split('/')[-1].replace('.faa', '') : path for path in Path(args.faas).rglob('*.faa')}
    fna_pathlists = {str(path).split('/')[-1].replace('.fna', '') : path for path in Path(args.fnas).rglob('*.fna')}

    pathlist = []
    for genome, faa_path in faa_pathlists.items():
        fna_path = fna_pathlists.get(genome, None)
        if fna_path:
            pathlist.append((faa_path, fna_path))
    
    output_json = str(args.output)
    workers = int(args.p)

    if workers is None:
        workers = multiprocessing.cpu_count() - 1

    logger.info("Measuring {} genomes with {} CPUs".format(len(pathlist), workers))
    filepath_gen = ((Path(faa_path), Path(fna_path)) for faa_path, fna_path in pathlist)
    with multiprocessing.Pool(workers) as p:
        pipeline_gen = p.map(_mapping_wrapper, filepath_gen)
        genomic_properties = dict(pipeline_gen)
    logger.info("Measured {} genomes".format(len(genomic_properties.keys())))
    
    json.dump(genomic_properties, open(output_json, 'w'))