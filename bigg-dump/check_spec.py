import json

SPEC_PATH = '/Users/zaking/repos/escher/map_spec.json'

def conforms_to_spec(obj):
    with open(SPEC_PATH, 'r') as f:
        spec = json.load(f)
    check_r(obj, spec)
        
def check_r(o, spec):
    if isinstance(spec, dict):
        if spec.keys()[0]=="*":
            for k, v in o.iteritems():
                try:
                    assert isinstance(k, str)
                except AssertionError:
                    raise AssertionError('Bad key: %s' % k)
                check_r(v, spec.values()[0])
        else:
            for k, v in spec.iteritems():
                try:
                    assert k in o
                except AssertionError:
                    print o
                    raise AssertionError('Missing key: %s' % k)
                check_r(o[k], spec[k])
    elif isinstance(spec, list):
        for x in o:
            check_r(x, spec[0])
    elif isinstance(spec, str):
        if spec=='String':
            the_type = str
        elif spec=="Float":
            the_type = float
        elif spec=="Integer":
            the_type = int
        elif spec=="Boolean":
            the_type = bool
        elif spec!="*":
            raise Exception("Bad spec string: %s" % spec)
        
        try:
            assert isinstance(o, the_type)
        except AssertionError:
            raise AssertionError('Bad type: %s should be %s' % (o, the_type))
