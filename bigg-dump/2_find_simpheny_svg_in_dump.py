from sys import argv
from os import listdir
from os.path import join
import re
import gzip
import json
from operator import itemgetter
import cPickle as pickle

from theseus import id_for_new_id_style, load_model

def main():
    try:
        search_dir = argv[1]
        model_name = argv[2]
    except:
        raise Exception('Not enough arguments.')
    try:
        in_file = argv[3]
    except:
        in_file = None

    if in_file:
        with open(in_file) as f: svg = f.read()
        metabolite_count = len(re.findall(r"fill:rgb\(255,160,128\)", svg))
        print 'metabolite_count: %d' % metabolite_count
        # ids = set(re.findall(r">([^<\[\]]*)(?:\[[a-z]+\])?</text>", svg))

    model = load_model(model_name)
    ids = [id_for_new_id_style(x.id) for x in model.reactions]
    
    scores = []
    for path in listdir(search_dir):
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
        # num_matches = len(re.findall(r"(%s)" % ("|".join(ids)), f.read()))
        # num_matches = 0
        scores.append((path, float(num_matches)/len(reactions), met_count, reaction_count))
        f.close()
    scores = sorted(scores, key=itemgetter(2), reverse=True)
    scores = sorted(scores, key=itemgetter(1))
    for score in scores:
        print score
    print ('path', 'id match %', 'node count', 'reaction_count')
    if in_file:
        print 'svg metabolite_count: %d' % metabolite_count    

if __name__=="__main__":
    main()
