#!/usr/bin/env python3

import argparse
from collections import defaultdict
import json
import logging
from pathlib import Path
import re

import numpy as np
import pandas as pd

class BacDiveData():
    """
    Parses BacDive API output to return reported information
    
    The BacDive API returns, for each BacDive id, a nested dictionary
    with information about the strain. The types of information can be
    found at https://api.bacdive.dsmz.de/strain_fields_information though
    the user should note that the '_' are replaced with ' ' and sections
    are keyed by the section name, not section ID. Multiple sets of values
    may be available for each type of data, in which case values are returned
    as a list of dictionaries instead of a dictionary.
    
    This class extracts information about strains: either info (taxid and
    genome accession) or conditions (pH, temperature, salinity, oxygen tolerance,
    media). Optimum values are the (average of) the reported optimum(s). Reported
    values refer to conditions in which positive growth was reported, which is for
    many strains a subset of the true range in which they can grow.
    """
    
    def __init__(self, entry):
        
        self.entry = entry

        self.strain_id = self.entry.get('General', {}).get( 'BacDive-ID', None)
        self.taxid_ncbi = self.get_taxid_ncbi()
        self.genome_accession_ncbi = self.get_genome_accession_ncbi()

        self.reported_media = self.get_reported_media()
        self.reported_temperatures = self.get_reported_temperatures()
        self.reported_phs = self.get_reported_phs()
        self.reported_salinities = self.get_reported_salinities()
        self.reported_oxygen_tolerances = self.get_reported_oxygen_tolerances()
        
        self.optimum_ph = self.get_optimum_ph()
        self.optimum_temperature = self.get_optimum_temperature()
        self.midpoint_salinity = self.compute_midpoint_salinity()
    
    def get_reported_media(self) -> list:
        subsection = self.entry.get('Culture and growth conditions', None).get('culture medium', {})
        media_ids = self._query_list_of_dicts(subsection, '@ref', 'growth', ['yes', 'positive'])
        return set(media_ids)
    
    def get_reported_temperatures(self) -> list:
        temperatures = []
        subsection = self.entry.get('Culture and growth conditions', {}).get('culture temp', {})
        for val in self._query_list_of_dicts(subsection, 'temperature', 'growth', ['yes', 'positive']):
            temperatures.extend(self._format_values(val))
        return set(temperatures)

    def get_reported_phs(self) -> list:
        phs = []
        subsection = self.entry.get('Culture and growth conditions', {}).get('culture pH', {})
        for val in self._query_list_of_dicts(subsection, 'pH', 'ability', ['yes', 'positive']):
            phs.extend(self._format_values(val))
        return set(phs)

    def get_reported_salinities(self) -> list:
        salinities = []
        subsection = self.entry.get('Physiology and metabolism', None).get('halophily', {})
        if isinstance(subsection, dict):
            salinities.extend(self.parse_halophily_dict(subsection))
        elif isinstance(subsection, list):
            for _dict in subsection:
                salinities.extend(self.parse_halophily_dict(_dict))
        return set(salinities)
    
    def get_optimum_ph(self) -> list:
        phs = []
        subsection = self.entry.get('Culture and growth conditions', {}).get('culture pH', {})
        for val in self._query_list_of_dicts(subsection, 'pH', 'type', ['optimum']):
            phs.extend(self._format_values(val))
        if len(phs) > 0:
            return np.mean(phs)
        else:
            return None
    
    def get_optimum_temperature(self) -> list:
        temperatures = []
        subsection = self.entry.get('Culture and growth conditions', {}).get('culture temp', {})
        for val in self._query_list_of_dicts(subsection, 'temperature',  'type', ['optimum']):
            temperatures.extend(self._format_values(val))
        if len(temperatures) > 0:
            return np.mean(temperatures)
        else:
            return None
        
    def get_reported_oxygen_tolerances(self) -> list:
        subsection = self.entry.get('Physiology and metabolism', None).get('oxygen tolerance', {})
        tolerances = self._query_list_of_dicts(subsection, 'oxygen tolerance', '', [None])
        return set(tolerances)
    
    def get_genome_accession_ncbi(self) -> str:
        subsection = self.entry.get('Sequence information', {}).get('Genome sequences', {})
        accessions = self._query_list_of_dicts(subsection, 'accession',  'database', ['ncbi'])
        if len(accessions) > 0:
            return accessions[0]
        else:
            return None
        
    def get_taxid_ncbi(self) -> str:
        """Returns taxid for lowest taxonomic level"""
        subsection = self.entry.get('General', {}).get('NCBI tax id', {})
        for level in ['strain', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'domain']:
            taxid = self._query_list_of_dicts(subsection, 'NCBI tax id',  'Matching level', [level])
            if taxid:
                return taxid[0]
        return None
    
    def compute_midpoint_salinity(self):
        if len(self.reported_salinities) > 0:
            return np.mean([min(self.reported_salinities), max(self.reported_salinities)])

    def parse_halophily_dict(self, halophily : dict):
        """
        Returns growth range of salinity (% NaCl) with the 
        following assumptions:
        
        - No growth > value means value is maximum, minimum is 0
        - Positive growth < value means value is maximum, minimum is 0
        - Other tested relations as point values for positive or inconsistent
            growth but not no growth
        """
        concentrations = []
        concentration = halophily.get('concentration', '')
        growth = halophily.get('growth')
        salt = halophily.get('salt')
        
        # conversions for NaCl only
        if 'g/L' in concentration:
            conversion = 0.1
        elif 'M' in concentration:
            conversion = 58.443 / 10.
        elif '%' in concentration:
            conversion = 1
        else:
            conversion = 1
            
        if salt == 'NaCl':
            if concentration.startswith('>') and growth == 'no':
                concentrations = [0.] + self._format_values(concentration)
            elif concentration.startswith('<') and growth in ['positive', 'inconsistent']:
                concentrations = [0.] + self._format_values(concentration)
            elif growth in ['positive', 'inconsistent']:
                concentrations = self._format_values(concentration)

        percent_nacl = [val * conversion for val in concentrations] 
        return [val for val in percent_nacl if val < 39]
    
    def _query_list_of_dicts(self, obj, key, required_key, required_vals : list):
        """
        Keys values from 1 dict or a list of dicts with a 
        condition that another key-value pair is present
        """
        arr = []
        if isinstance(obj, dict):
            val = self._query_dict_conditionally(obj, key, required_key, required_vals)
            if val:
                arr.append(val)
        elif isinstance(obj, list): # multiple values
            for subobj in obj:
                if isinstance(subobj, dict):
                    val = self._query_dict_conditionally(subobj, key, required_key, required_vals)
                    if val:
                        arr.append(val)
        return arr
    
    def _query_dict_conditionally(self, obj : dict, key, required_key, required_vals : list):
        """Lookup value if other value in dictionary among accepted values"""
        if obj.get(required_key, None) in required_vals:
            return obj.get(key, None)
        else:
            return None

    def _format_values(self, string : str) -> list:
        """
        Uses regex and replace to extract non-float characters from
        strings and correct for typos in data entry. If a range,
        like 3.4-8.4, both values will be returned. Otherwise, one
        value will be returned.
        """
        regex = re.compile(r"[-+]?(?:\d*\.*\d+)")
        if '-' in string:
            return [float(regex.search(val).group(0).replace('..', '.')) for val in string.split('-') if len(val) > 0]
        else:
            return [float(regex.search(string).group(0))]
        



class MakeFeatureTable():
    """
    Creates feature table from genomic data and trait data.
    
    Memory intensive due to the size of the data and loading
    the JSONs into memory at once.
    """
    
    def __init__(self, trait_data : Path,  genomic_data : Path,
                 genomic_metadata : list, tsv_path : Path,
                ):
        
        self.tsv_path = tsv_path
        self.genomic_metadata_paths = genomic_metadata
        with open(trait_data) as fh:
            self.trait_data = json.loads(fh.read())
        with open(genomic_data) as fh:
            self.genomic_data = json.loads(fh.read())
        
        self.genomic_metadata = None
        self.genomic_features = None
        self.trait_features = None
        self.features = None
    
    def create_genomic_features(self, data : dict):
        """
        Loads genomic data to features. 
        """
        features = {}
        
        for protein_set in ['all', 'extracellular_soluble']:            
            for k, v in data.items():
                if k.startswith('histogram_'):
                    # Convert counts to frequency
                    total_counts = sum(v.values())
                    for subk, subv in v.items():
                        features[protein_set + '_' + k.replace('histogram_', '') + '_' + subk] = float(subv / total_counts)
                else:
                    features[protein_set + '_' + k] = v
        return features
            
    def create_trait_features(self, data : dict, source : str='bacdive'):
        """
        Loads trait data, differently by source, to a set of
        features that will be used for modeling.
        
        Accepted sources are currently only 'bacdive'
        """
        if source == 'bacdive':
            strain = BacDiveData(data)
        else:
            raise ValueError("Source must be one of: ['bacdive']")

        features = {
            'ncbi_accession' : strain.genome_accession_ncbi,
            'ncbi_taxid' : strain.taxid_ncbi,
            'strain_id' : strain.strain_id,
            'optimum_ph' : strain.optimum_ph,
            'optimum_temperature' : strain.optimum_temperature,
            'midpoint_salinity' : strain.midpoint_salinity,
                   }

        features.update(self.onehot_range(strain.reported_salinities, 0, 38.4, 0.5, 'nacl_')) # salinity range
        features.update(self.onehot_range(strain.reported_phs, 0, 14, 0.25, 'ph_')) # pH range 
        features.update(self.onehot_range(strain.reported_temperatures, 0, 100, 1, 'temp_')) # temperature range
        features.update(self.onehot_oxygen_tolerance(strain.reported_oxygen_tolerances)) # o2 tolerance

        return features
    
    def onehot_range(self, arr, min_bin : float, max_bin : float, step : float, prefix : str) -> dict:
        onehot_dict = {}
        for bin_floor in np.arange(min_bin, max_bin, step):
            feature_name = prefix + str(bin_floor)
            if min(arr, default=-1000) <= bin_floor <= max(arr, default=-1000):
                onehot_dict[feature_name] = 1
            else:
                onehot_dict[feature_name] = 0
        return onehot_dict
    
    def onehot_oxygen_tolerance(self, tolerances : set) -> dict:
        """
        Returns a dictionary of oxygen tolerance definitions with 1
        indicating the organism has that tolerance. Subtypes of aerobe 
        and anaerobe lead to assignment of 1 to aerobe or anaerobe, 
        respectivey, with facultative anaerobes only assigned to aerobe. 
        """
        onehot_tolerances = {'aerobe' : None, 
                       'anaerobe' : None, 
                       'microaerophile' : None , 
                       'facultative anaerobe' : None, 
                       'obligate aerobe' : None, 
                       'obligate anaerobe' : None, 
                       'facultative aerobe' : None, 
                       'aerotolerant' : None, 
                       'microaerotolerant' : None}
        
        aerobe_subtypes = {'facultative anaerobe', 'obligate aerobe', 'facultative aerobe', 'microaerophile'}
        anaerobe_subtypes = {'obligate anaerobe', 'facultative aerobe'}
        
        for tolerance in tolerances:
            onehot_tolerances[tolerance] = 1
        if len(tolerances.intersection(aerobe_subtypes)) > 0:
            onehot_tolerances['aerobe'] = 1
        if len(tolerances.intersection(anaerobe_subtypes)) > 0:
            onehot_tolerances['anaerobe'] = 1
        
        return onehot_tolerances
        

    def load_genomic_metadata(self, source : str='gtdb'):
        """
        Loads and formats GTDB metadata. Only essential colums are kept.
        Column 'ncbi_acccession' and 'ncbi_taxid' matches information stored
        by BacDive.
        """
        
        if self.genomic_metadata is None:
            if source == 'gtdb':
                sep = '\t'
                retain_columns = {'accession' : str, 'ncbi_taxid' : int, 'ncbi_genbank_assembly_accession' : str,
                                  'ncbi_taxonomy' : str,  'gtdb_taxonomy' : str,  'gtdb_genome_representative' : str, 
                                  'coding_density' : float, 'gc_percentage' : float, 'genome_size' : int, 
                                 }

                df = pd.concat([pd.read_csv(tsv, sep=sep, header=0, dtype=retain_columns) for tsv in self.genomic_metadata_paths], axis=0).loc[:, retain_columns.keys()]

                # Match BacDive format
                df['ncbi_accession'] = df['ncbi_genbank_assembly_accession'].str.split('.', expand=True)[0].tolist()
                df = df.drop(columns=['ncbi_genbank_assembly_accession'])
                df = df.set_index('accession')
                self.genomic_metadata = df.sort_values('genome_size', ascending=False)
            else:
                raise ValueError("Source must be one of: ['gtdb']")
            
        return self.genomic_metadata

    def load_features(self):
        """
        Loads a dictionary of features for trait data, genomic data,
        and metadata. To be memory efficient, one strain is loaded at
        a time. Each strain is assigned to the genome of the GTDB
        representative or discarded. A representative genome can be assigned
        to multiple strains. 
        """
        self.genomic_metadata = self.load_genomic_metadata(source='gtdb')
        
        # strain_accessions = set()
        # for strain_id, data in self.trait_data.items():
        #     trait_features = self.create_trait_features(data, source='bacdive')
        #     strain_accessions.add(trait_features['ncbi_accession'])

        accession_dict = self.genomic_metadata.reset_index().set_index('ncbi_accession')['gtdb_genome_representative'].to_dict()
        
        if self.features == None:
            self.features = {}
            n = 0
            for strain_id, data in self.trait_data.items():
                strain_features = {}
                trait_features = self.create_trait_features(data, source='bacdive')
                representative_key = accession_dict.get(trait_features['ncbi_accession'], None)
                if representative_key:
                    representative_metadata = self.genomic_metadata.loc[[representative_key]].head(1).to_dict(orient='records')[0] # first row if dup

                    strain_features.update(representative_metadata)
                    strain_features.update(self.create_genomic_features(self.genomic_data[representative_key]))
                    strain_features.update(trait_features)
                    self.features[strain_id] = strain_features
                        
        return self.features
    
    def write_feature_table(self):
        pd.DataFrame.from_dict(self.load_features(), orient='index').to_csv(self.tsv_path, sep='\t', index=None)
        return self.tsv_path



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
                prog='',
                description=''
                )
    
    parser.add_argument('-b', '--bacdive-json', help='JSON output of query_bacdive', required=True)

    args = parser.parse_args()
    
    print(extract_features_from_bacdive_data(data_json=args.bacdive_json))