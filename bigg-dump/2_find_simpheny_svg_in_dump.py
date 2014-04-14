from sys import argv
from os import listdir
from os.path import join
import re
import gzip
import json
from operator import itemgetter

def main():
    try:
        in_file = argv[1]
        search_dir = argv[2]
    except:
        raise Exception('Not enough arguments.')

    with open(in_file) as f:
        svg = f.read()

    metabolite_count = len(re.findall(r"fill:rgb\(255,160,128\)", svg))
    print metabolite_count
    ids = set(re.findall(r">([^<\[\]]*)(?:\[[a-z]+\])?</text>", svg))

    scores = []
    for path in listdir(search_dir):
        if path.endswith('.gz'):
            f = gzip.open(join(search_dir, path), 'r')
        else:
            f = open(join(search_dir, path), 'r')
        # compare the metabolite count 
        m = json.load(f)
        try:
            diff = abs(len(m['MAPNODE']) - metabolite_count)
        except KeyError:
            continue
        # compare the id matches
        f.seek(0)
        num_matches = len(re.findall(r"(%s)" % ("|".join(ids)), f.read()))
        scores.append((path, num_matches, diff))
        f.close()
    scores = sorted(scores, key=itemgetter(1))
    scores = sorted(scores, key=itemgetter(2), reverse=True)
    for score in scores:
        print score
    print ('path', 'id matches', 'difference in metabolite count')

if __name__=="__main__":
    main()
