"""Process svg files to generate escher json files.

"""
import sys
import gzip
import json
from pprint import pprint
from math import atan2, isnan
from itertools import izip
from numpy import inf
from os.path import join, basename
from os import listdir
import logging
import pandas as pd
import cobra
import cPickle as pickle
import jsonschema
from urllib2 import urlopen

from check_spec import conforms_to_spec
from theseus import id_for_new_id_style, load_model

def main():
    """Load a processed dump file (.json or .json.gz), and generate a map for the
    builder tool.

    """
    try:
        in_files = sys.argv[1:-2]
        model_name = sys.argv[-2]
        out_directory = sys.argv[-1]
    except IndexError:
        raise Exception("Not enough arguments")

    if len(in_files)==0:
        raise Exception("Not enough arguments")

    if in_files[0].endswith('.tsv'):
        in_files = (pd.DataFrame
                    .from_csv(in_files[0], sep='\t')
                    .loc[:, 'path'])

    size = len(in_files)
    for i, filename in enumerate(in_files):
        # progress
        sys.stdout.write('\r')
        # the exact output you're looking for:
        sys.stdout.write("%d / %d" % (i + 1, size))
        sys.stdout.flush()
        
        save_map(filename, out_directory, model_name)
            
def save_map(filename, out_directory, model_name):
    
    if filename.endswith('.json.gz'):
        with gzip.open(filename, "r") as f:
            data = json.load(f)
        out_file = join(out_directory, basename(filename).replace('.json.gz', '_map.json'))
    elif filename.endswith('.json'):
        with open(filename, "r") as f:
            data = json.load(f)
        out_file = join(out_directory, basename(filename).replace('.json', '_map.json'))
    else:
        logging.warn('Not loading file %s' % filename)

    # get the cobra model
    model = load_model(model_name)

    # get the compartment dictionary
    df = pd.DataFrame.from_csv("compartment_id_key.csv")
    compartment_id_dictionary = {}
    for row in df.itertuples(index=True):
        compartment_id_dictionary[row[0]] = row[1:3]
        
    # major categories
    reactions = []; line_segments = []; text_labels = []; nodes = []
    for k, v in data.iteritems():
        if k=="MAPREACTION": reactions = v
        elif k=="MAPLINESEGMENT": segments = v
        elif k=="MAPTEXT": text_labels = v
        elif k=="MAPNODE": nodes = v
        else: raise Exception('Unrecognized category: %s' % k)
        
    # do the nodes
    nodes = parse_nodes(nodes, compartment_id_dictionary)

    # do the segments
    parse_segments(segments, reactions, nodes)
        
    # do the reactions
    reactions = parse_reactions(reactions, model, nodes)

    # do the text labels
    text_labels = parse_labels(text_labels)

    # compile the data
    out = {}
    out['nodes'] = nodes
    out['reactions'] = reactions
    out['text_labels'] = text_labels

    # translate everything so x > 0 and y > 0
    # out = translate_everything(out)
    
    # for export, only keep the necessary stuff
    node_keys_to_keep = ['node_type', 'x', 'y', 'name', 'bigg_id', 'label_x',
                         'label_y', 'node_is_primary', 'connected_segments']
    segment_keys_to_keep = ['from_node_id', 'to_node_id', 'b1', 'b2']
    reaction_keys_to_keep = ['segments', 'name', 'reversibility',
                             'bigg_id', 'label_x', 'label_y', 'metabolites',
                             'gene_reaction_rule']
    text_label_keys_to_keep = ['x', 'y', 'text']    
    for k, node in out['nodes'].iteritems():
            only_keep_keys(node, node_keys_to_keep)
    for k, reaction in out['reactions'].iteritems():
        if 'segments' not in reaction: continue
        for k, segment in reaction['segments'].iteritems():
            only_keep_keys(segment, segment_keys_to_keep)
        only_keep_keys(reaction, reaction_keys_to_keep)
    for k, text_label in out['text_labels'].iteritems():
        only_keep_keys(text_label, text_label_keys_to_keep)

    # get max width and height
    min_max = {'x': [inf, -inf], 'y': [inf, -inf]}
    for node in nodes.itervalues():
        if node['x'] < min_max['x'][0]: min_max['x'][0] = node['x']
        if node['x'] > min_max['x'][1]: min_max['x'][1] = node['x']
        if node['y'] < min_max['y'][0]: min_max['y'][0] = node['y']
        if node['y'] > min_max['y'][1]: min_max['y'][1] = node['y']
    width = min_max['x'][1] - min_max['x'][0]
    height = min_max['y'][1] - min_max['y'][0]
    out['canvas'] = { 'x': min_max['x'][0] - 0.05 * width,
                      'y': min_max['y'][0] - 0.05 * height,
                      'width': width + 0.10 * width,
                      'height': height + 0.10 * height}

    header = {
        "schema": "https://zakandrewking.github.io/escher/escher/jsonschema/1-0-0#",
        "homepage": "https://zakandrewking.github.io/escher",
        "map_id": basename(filename).replace('.json', '').replace('.gz', ''),
        "map_name": "",
        "map_description": ""
        }

    the_map = [header, out]
        
    # make sure it conforms
    # schema = json.loads(urlopen('https://zakandrewking.github.io/escher/schema/v1/schema')
    #                     .read())
    with open('/home/king/repos/escher/escher/jsonschema/1-0-0', 'r') as f:
        schema = json.load(f)
    jsonschema.validate(the_map, schema)
    print 'Map is valid'
    
    with open(out_file, 'w') as f: json.dump(the_map, f, allow_nan=False)

def parse_nodes(nodes, compartment_id_key):
    for node in nodes:
        # assign new keys
        try_assignment(node, 'MAPOBJECT_ID', 'object_id',
                       cast=str, require=True)
        try_assignment(node, 'MAPNODENODETYPE', 'node_type',
                       cast=lambda x: str(x).lower(), require=True)
        try_assignment(node, 'MAPNODEPOSITIONX', 'x',
                       cast=float, require=True)
        try_assignment(node, 'MAPNODEPOSITIONY', 'y',
                       cast=float, require=True)
        if node['node_type'] == 'metabolite':
            try_assignment(node, 'MAPNODELABELPOSITIONX', 'label_x',
                           cast=float, fallback=node['x'])
            try_assignment(node, 'MAPNODELABELPOSITIONY', 'label_y',
                           cast=float, fallback=node['y'])
            try_assignment(node, 'MOLECULEABBREVIATION', 'bigg_id',
                           cast=lambda x: id_for_new_id_style(x, is_metabolite=True),
                           fallback='')
            try_assignment(node, 'MOLECULEOFFICIALNAME', 'name',
                           cast=str, fallback='')
            try_assignment(node, 'MAPNODEISPRIMARY', 'node_is_primary',
                           cast=lambda x: True if x=='Y' else False, fallback=False)
            try_assignment(node, 'MAPNODECOMPARTMENT_ID', 'compartment_id', cast=int)
            try_assignment(node, 'MAPNODECOMPARTMENT_ID', 'compartment_name',
                           cast=lambda x: compartment_id_key[int(x)][0])
            try_assignment(node, 'MAPNODECOMPARTMENT_ID', 'compartment_letter',
                           cast=lambda x: compartment_id_key[int(x)][1])
            if node['bigg_id'] != '' and node['compartment_letter'] is not None:
                node['bigg_id'] = "%s_%s" % (node['bigg_id'],
                                             node['compartment_letter'])
        
    # Make into dictionary
    return {a['object_id']: a for a in nodes
            if a['node_type'] in ['metabolite', 'multimarker', 'midmarker']}

def check_and_add_to_nodes(nodes, node_id, segment_id, reaction_id):
    try:
        node = nodes[node_id]
    except KeyError:
        return None
    if 'connected_segments' not in node:
        node['connected_segments'] = []
    node['connected_segments'].append({'reaction_id': reaction_id,
                                       'segment_id': segment_id})
    return node_id
        
def to_x_y(array):
    x = float(array[0])
    y = float(array[1])
    if isnan(x) or isnan(y):
        return None
    return {'x': x, 'y': y}

def parse_segments(segments, reactions, nodes):
    # do the segments
    segment_id = 0
    for segment in segments:
        # get the reaction
        reaction_id = str(segment['MAPLINESEGMENTREACTION_ID'])
        reaction = [a for a in reactions if str(a['MAPOBJECT_ID'])==reaction_id]

        if len(reaction) > 1: reaction = reaction[0] # raise Exception('Too many matches')
        else: reaction = reaction[0]
            
        # get the nodes, and make sure they both exist before adding
        from_node_id = check_and_add_to_nodes(nodes, segment['MAPLINESEGMENTFROMNODE_ID'],
                                              str(segment_id), str(reaction['MAPOBJECT_ID']))
        to_node_id = check_and_add_to_nodes(nodes, segment['MAPLINESEGMENTTONODE_ID'],
                                            str(segment_id), str(reaction['MAPOBJECT_ID']))
        if from_node_id is not None and to_node_id is not None:
            segment['from_node_id'] = from_node_id
            segment['to_node_id'] = to_node_id
            try:
                segment['b1'] = to_x_y(segment['MAPLINESEGMENTCONTROLPOINTS'][0])
                segment['b2'] = to_x_y(segment['MAPLINESEGMENTCONTROLPOINTS'][1])
            except KeyError:
                segment['b1'] = None
                segment['b2'] = None

            if 'segments' not in reaction:
                reaction['segments'] = {}
            reaction['segments'][str(segment_id)] = segment
            segment_id = segment_id + 1

def parse_reactions(reactions, model, nodes):
    d = {'x': 10, 'y': -10}
    for reaction in reactions:
        try_assignment(reaction, 'REACTIONOFFICIALNAME', 'name',
                       cast=str, require=True)
        try_assignment(reaction, 'REACTIONABBREVATION', 'bigg_id',
                       cast=id_for_new_id_style, require=True)
        try_assignment(reaction, 'MAPOBJECT_ID', 'object_id',
                       cast=str, require=True)
            
        # get reaction label position. set to the midmarker
        try:
            connected_nodes = set(reduce(lambda x, y: x + [y['to_node_id'], y['from_node_id']],
                                         reaction['segments'].itervalues(), []))
        except KeyError:
            continue
        
        for n in connected_nodes:
            node = nodes[n]
            if node['node_type']=='midmarker':
                reaction['label_x'] = node['x'] + d['x']
                reaction['label_y'] = node['y'] + d['y']
        
        # use the cobra model
        try:
            cobra_reaction = model.reactions.get_by_id(reaction['bigg_id'])
            if (cobra_reaction.lower_bound < 0 and cobra_reaction.upper_bound <= 0):
                # reverse the reaction
                reaction['reversibility'] = False
                reaction['metabolites'] = [{'bigg_id': k.id, 'coefficient': -v}
                                           for k, v in cobra_reaction._metabolites.iteritems()]
            else:
                # get the reversibility
                reaction['reversibility']  = (cobra_reaction.lower_bound < 0)
                # get the metabolites
                reaction['metabolites'] = [{'bigg_id': k.id, 'coefficient': v}
                                            for k, v in cobra_reaction._metabolites.iteritems()]
            reaction['gene_reaction_rule'] = cobra_reaction.gene_reaction_rule
        except KeyError:
            print 'Could not add reaction %s' % reaction['bigg_id']
            del reaction['segments']
        
    # take out the reactions without segments
    return {r['object_id']: r for r in reactions if 'segments' in r}

def parse_labels(labels):
    for label in labels:
        try_assignment(label, 'MAPTEXTPOSITIONX', 'x',
                       cast=float, require=True)
        try_assignment(label, 'MAPTEXTPOSITIONY', 'y',
                       cast=float, require=True)
        try_assignment(label, 'MAPOBJECT_ID', 'object_id',
                       cast=str, require=True)
        try_assignment(label, "MAPTEXTCONTENT", 'text',
                       cast=str, require=True)
    return {r['object_id']: r for r in labels}

def only_keep_keys(d, keys):
    for k, v in d.items():
        if k not in keys:
            del d[k]

def try_assignment(node, key, new_key, cast=None, fallback=None, require=False):
    try:
        if cast is None:
            node[new_key] = node[key]
        else:
            node[new_key] = cast(node[key])
        try:
            assert (not isnan(node[new_key]))
        except TypeError:
            pass
    except (KeyError, AssertionError) as err:
        if require:
            raise err
        else:
            node[new_key] = fallback

def translate_everything(out):
    def get_min(a):
        def check(d, mins):
            if type(d) is not dict: return
            
            xs = []; ys = []
            if 'x' in d and d['x'] is not None:
                xs.append(d['x'])
            if 'label_x' in d and d['label_x'] is not None:
                xs.append(d['label_x'])
            if 'y' in d and d['y'] is not None:
                ys.append(d['y'])
            if 'label_y' in d and d['label_y'] is not None:
                ys.append(d['label_y'])

            mins['x'] = min(xs + [mins['x']])
            mins['y'] = min(ys + [mins['y']])

            for k, v in d.iteritems():
                check(v, mins)

        mins = {'x': 0, 'y': 0}
        check(a, mins)
        return mins['x'], mins['y']
    
    def translate(d, subtract_x, subtract_y):
        if type(d) is not dict: return d

        if 'x' in d and d['x'] is not None:
            d['x'] = d['x'] - subtract_x
        if 'label_x' in d and d['label_x'] is not None:
            d['label_x'] = d['label_x'] - subtract_x
        if 'y' in d and d['y'] is not None:
            d['y'] = d['y'] - subtract_y
        if 'label_y' in d and d['label_y'] is not None:
            d['label_y'] = d['label_y'] - subtract_y
        
        for k, v in d.iteritems():
            d[k] = translate(v, subtract_x, subtract_y)
            
        return d

    subtract_x, subtract_y = get_min(out)
    print subtract_x, subtract_y
    return translate(out, subtract_x, subtract_y)

if __name__=="__main__":
    main()
