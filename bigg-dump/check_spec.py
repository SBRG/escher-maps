import json

SPEC_PATH = '/Users/zaking/www/escher-private/escher/map_spec.json'

def conforms_to_spec(obj):
    with open(SPEC_PATH, 'r') as f:
        spec = json.load(f)
    check_r(obj, spec["spec"], spec["can_be_none"])
        
def check_r(o, spec, can_be_none):
    if isinstance(spec, dict):
        if spec.keys()[0]=="*":
            for k, v in o.iteritems():
                if v is None and k in can_be_none:
                    continue
                try:
                    assert isinstance(k, str)
                except AssertionError:
                    raise AssertionError('Bad key: %s' % k)
                check_r(v, spec.values()[0], can_be_none)
        else:
            for k, v in spec.iteritems():
                try:
                    assert k in o
                except AssertionError:
                    print o
                    raise AssertionError('Missing key: %s' % k)
                if o[k] is None and k in can_be_none:
                    continue
                check_r(o[k], spec[k], can_be_none)
    elif isinstance(spec, list):
        for x in o:
            check_r(x, spec[0], can_be_none)
    elif isinstance(spec, str) or isinstance(spec, unicode):
        if spec=='String':
            the_type = (str, unicode)
        elif spec=="Float":
            the_type = (float, int)
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
    else:
        raise Exception('Bad spec value: %s with type %s' % (spec, type(spec)))
