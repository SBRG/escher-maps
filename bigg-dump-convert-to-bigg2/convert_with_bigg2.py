#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Upgrade and convert maps to match BiGG 2.

Convert the BiGG 1 dumps to valid Escher maps, and assign IDs from the BiGG 2
database.

NOTE: This requires Escher 1.2-dev.

"""

# locations
map_paths = ['maps', '../1-0-0/maps']
model_id_mapping = {'EcoliCore': 'e_coli_core',
                    'Recon1': 'RECON1'}


import requests
from os import makedirs, listdir
from os.path import join, isdir
import json
import cobra
import re
import logging
import sys

from escher.convert_map import convert
import ome
from ome.models import *
from ome.base import *


# configure logger
logging.basicConfig(stream=sys.stdout, level=logging.INFO)


class NotFoundError(Exception):
    pass

def make_url(path):
    """make a bigg api request url"""
    return 'http://bigg.ucsd.edu/api/v2/%s' % path.lstrip('/')

def parse_model(filename):
    """Get the model from the filename."""
    if '_' in filename:
        return filename.split('_', 1)[0]
    else:
        return filename.split('.', 1)[0]

def parse_map_name(filename):
    """Get the map name from the filename."""
    return filename.replace('.json', '')

def load_bigg_model(model_id):
    """Try to download this model from BiGG."""
    # check id mapping
    try:
        model_id = model_id_mapping[model_id]
    except KeyError:
        pass
    
    # download the model
    url = make_url('/models/%s/download' % model_id)
    response = requests.get(url)
    if response.status_code != 200:
        raise NotFoundError('Could not find %s (url %s, code %d)' %
                            (model_id, url, response.status_code))
    return cobra.io.json._from_dict(response.json())

def fix_filename(filename):
    """Make an Escher compatible filename."""
    return '.'.join(filename.split('_', 1))

def fix_mapping(mapping, mapping_type):
    """Some published models to not match maps. Apply this to a mapping dictionary
    to include specific fixes.

    """
    subs = [
        (r'^([DL])_', r'\1-'), # in published model, in map
        (r'_([CDLR0-9])(_|$)', r'-\1\2'),
        (r'_DASH_', r'-'),
        (r'_LPAREN_(\w+)_RPAREN_', r'(\1)'),
        (r"'", ''),
        (r'_([a-z])_?$', r'(\1)'),
        (r'__atp_', r'(atp)'),
        (r'__gtp_', r'(gtp)'),
        (r'_gtp_', r'(gtp)'),
        (r'__adp_', r' (adp)'),
        (r'_nadp_', r'(nadp)'),
        (r'_10$', r'_1.0'),
        (r'_FMN', r' (FMN)'),
    ]

    for regex, replace in subs:
        for k, v in list(mapping.items()):
            modified_old_key = re.sub(regex, replace, k)
            if modified_old_key != k:
                mapping[modified_old_key] = v

    return mapping

def get_reaction_mapping(model_id):
    """Get the old reaction to new reaction mapping as a dictionary (old are keys)."""
    session = Session()
    result_db = (session.query(Synonym.synonym, Reaction.bigg_id)
                 .join(OldIDSynonym, OldIDSynonym.synonym_id == Synonym.id)
                 # .filter(OldIDSynonym.type == 'model_reaction') # for debugging
                 .join(ModelReaction, ModelReaction.id == OldIDSynonym.ome_id)
                 .join(Reaction, Reaction.id == ModelReaction.reaction_id)
                 .join(Model, Model.id == ModelReaction.model_id)
                 .filter(Model.bigg_id == model_id)
                 .all())

    # find duplicates
    seen = set(); duplicates = set(); mapping = {}
    for _, new in result_db:
        if new not in seen:
            seen.add(new)
        else:
            duplicates.add(new)

    for old, new in result_db:
        # for duplicates, make them all _copy1
        if new in duplicates:
            mapping[old] = '%s_copy1' % new
        else:
            mapping[old] = new
        
    session.close()
    return fix_mapping(mapping, 'reaction')


def get_metabolite_mapping(model_id):
    """Get the old metabolite to new metabolite mapping as a dictionary (old are keys)."""
    session = Session()
    result_db = (session.query(Synonym.synonym, Component.bigg_id, Compartment.bigg_id)
                 .join(OldIDSynonym, OldIDSynonym.synonym_id == Synonym.id)
                 # for debugging
                 # .filter(OldIDSynonym.type == 'model_compartmentalized_component') 
                 .join(ModelCompartmentalizedComponent,
                       ModelCompartmentalizedComponent.id == OldIDSynonym.ome_id)
                 .join(CompartmentalizedComponent,
                       CompartmentalizedComponent.id == ModelCompartmentalizedComponent.compartmentalized_component_id)
                 .join(Component,
                       Component.id == CompartmentalizedComponent.component_id)
                 .join(Compartment,
                       Compartment.id == CompartmentalizedComponent.compartment_id)
                 .join(Model,
                       Model.id == ModelCompartmentalizedComponent.model_id)
                 .filter(Model.bigg_id == model_id)
                 .all())
    mapping = {x[0]: '%s_%s' % x[1:3] for x in result_db}
    session.close()
    return fix_mapping(mapping, 'metabolite')


def get_dirs(path_list):
    files = []; outputs = set()
    for path in path_list:
        # recurse
        for filename in listdir(path):
            filepath = join(path, filename) 
            if isdir(filepath):
                files = files + get_dirs([filepath])
            else:
                output = '%s-converted' % path
                outputs.add(output)
                files.append((filepath, filename, output))
    # make the directories
    for output in outputs:
        try:
            makedirs(output)
        except OSError:
            pass
    return files


def main():
    # go though the maps
    for filepath, filename, output_path in get_dirs(map_paths):
        logging.info('Converting %s' % filename)
        # # debug
        # if 'iAF692' not in filename:
        #     continue
        
        # get the model
        model_id = parse_model(filename)
        try:
            model = load_bigg_model(model_id)
        except NotFoundError as e:
            print(filepath, e)
            continue

        # get the id mapping dictionaries. Eventually, this should be available
        # in the API.
        reaction_id_mapping = get_reaction_mapping(model.id)
        metabolite_id_mapping = get_metabolite_mapping(model.id)
        gene_id_mapping = None # get_gene_mapping(model_id)

        # convert the map
        with open(filepath, 'r') as f:
            map_in = json.load(f)
        map_out = convert(map_in, model, map_name=parse_map_name(filename),
                          reaction_id_mapping=reaction_id_mapping,
                          metabolite_id_mapping=metabolite_id_mapping,
                          gene_id_mapping=gene_id_mapping)

        # write
        with open(join(output_path, fix_filename(filename)), 'w') as f:
            json.dump(map_out, f)

if __name__ == '__main__':
    main()
