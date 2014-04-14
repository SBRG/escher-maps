import sys
import json

def main():
    try:
        in_file = sys.argv[1]
    except IndexError:
        raise Exception("Not enough arguments")

    with open(in_file, 'r') as f:
        m = json.load(f)

    count = 0
    d = {'x': 10, 'y': -10}
    
    for k, v in m['reactions'].iteritems():
        connected_nodes = set(reduce(lambda x,y: x + [y['to_node_id'], y['from_node_id']],
                                     v['segments'].itervalues(), []))
        for n in connected_nodes:
            node = m['nodes'][n]
            if node['node_type']=='midmarker':
                m['reactions'][k]['label_x'] = node['x'] + d['x']
                m['reactions'][k]['label_y'] = node['y'] + d['y']
                count = count+1

    print 'Fixed %d labels' % count
                
    with open(in_file.replace('.json', '_realigned_labels.json'), 'w') as f:
        json.dump(m, f)

if __name__=="__main__":
    main()
