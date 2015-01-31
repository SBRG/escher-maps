import cobra.io
import json
from os import listdir
from os.path import join
from escher.convert_map import convert
import re

if __name__=="__main__":
    from sys import argv, exit
    if len(argv) < 2:
        print 'Usage: python update_gene_ids_in_maps.py model_id'
        exit()
        
    # find the model and the organism
    model_id = argv[1]
    model_dir = 'models'
    the_organism = None 
    for orgdir in listdir(model_dir):
        if the_organism is not None:
            break
        if orgdir.startswith('.'):
            continue
        orgpath = join(model_dir, orgdir)
        for model_file in listdir(orgpath):
            if model_file.startswith('.'):
                continue
            if model_file.split('.')[0] == model_id:
                the_organism = orgdir
                model_path = join(orgpath, model_file)
                model = cobra.io.load_json_model(model_path)
                gene_names = {gene.id.split('.')[0]: gene.name for gene in model.genes}
                # remember the gene names
                for r in model.reactions:
                    r.gene_reaction_rule = re.sub('\.[0-9]+', '', r.gene_reaction_rule)
                # remove old genes
                for old_id in gene_names.values():
                    model.genes.get_by_id(old_id).remove_from_model()
                # update the gene names with the old names
                for g in model.genes:
                    g.name = g.id
                model.id = model_id
                cobra.io.save_json_model(model, model_path)
                print 'Found model %s %s' % (the_organism, model_path)
                break
                
    if the_organism is None: 
        print 'Could not find model'
        exit()
                
    # apply model to maps
    map_org = join('maps', the_organism)
    for map_file in listdir(map_org):
        if map_file.startswith('.'):
            continue
        map_path = join(map_org, map_file)
        print 'Writing map %s' % map_path
        with open(map_path, 'r') as f:
            map = json.load(f)
        for r in map[1]['reactions'].itervalues():
            r['gene_reaction_rule'] = re.sub('\.[0-9]+', '', r['gene_reaction_rule'])
            for g in r['genes']:
                g['bigg_id'] = g['bigg_id'].split('.')[0]
        new_map = convert(map, model) # TODO the gene names are not being fixed here
        # import ipdb; ipdb.set_trace()
        with open(map_path, 'w') as f:
            json.dump(new_map, f)
