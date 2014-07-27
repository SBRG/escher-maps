"""Find map that matches the given model.

Returns a list of maps

"""
from sys import argv
from os import listdir
from os.path import join
import re
import gzip
import json
from operator import itemgetter
import cPickle as pickle
from time import sleep
import sys
import pandas as pd

from theseus import id_for_new_id_style, load_model

def main():
    try:
        search_dir = argv[1]
        model_name = argv[2]
    except:
        raise Exception('Not enough arguments.')

    model = load_model(model_name)
    ids = [id_for_new_id_style(x.id) for x in model.reactions]
    
    scores = []; size=len(listdir(search_dir))
    for i, path in enumerate(listdir(search_dir)):

        # progress
        sys.stdout.write('\r')
        # the exact output you're looking for:
        sys.stdout.write("%d / %d" % (i + 1, size))
        sys.stdout.flush()
                        
        if path.endswith('.gz'):
            f = gzip.open(join(search_dir, path), 'r')
        else:
            f = open(join(search_dir, path), 'r')
        # (1) Compare the metabolite count 
        m = json.load(f)
        try:
            met_count = len(m['MAPNODE'])
            reaction_count = len(m['MAPREACTION'])
            # diff = abs(len(m['MAPNODE']) - metabolite_count)
        except KeyError:
            continue
        # (2) Compare the reaction ids to the cobra model
        # f.seek(0)
        num_matches = 0
        try:
            reactions = m['MAPREACTION']
        except KeyError:
            continue
        for reaction in reactions:
            try:
                an_id = reaction['REACTIONABBREVATION'] 
            except KeyError:
                continue
            if id_for_new_id_style(an_id) in ids:
                num_matches = num_matches + 1
        # quit if not 100%
        if num_matches < len(reactions): continue
        scores.append((join(search_dir, path),
                       float(num_matches) / len(reactions),
                       met_count, reaction_count))
        f.close()
    scores = sorted(scores, key=itemgetter(2), reverse=True)
    scores = sorted(scores, key=itemgetter(1))
    outfile = '%s_maps.tsv' % model_name
    print
    print 'saving to %s' % outfile
    (pd
     .DataFrame(scores, columns=['path', 'score', 'n_metabolites', 'n_reactions'])
     .to_csv(outfile, sep='\t'))

if __name__=="__main__":
    main()
