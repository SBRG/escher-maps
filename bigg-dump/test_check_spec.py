from check_spec import conforms_to_spec
from pytest import raises

def test_conforms_to_spec():
    o = {"reactions": {"0": {"bigg_id": "GLCtex",
                             "reversibility": True,
                             "metabolites":
                             {"glc__D_p": {"coefficient":1,
                                           "index":0,
                                           "is_primary":True },
                              "glc__D_e": {"coefficient":-1,
                                           "index":0,
                                           "is_primary":True }
                          },
                             "label_x":3312.17,
                             "label_y":1910.74,
                             "name":"glucose transport via diffusion (extracellular to periplasm)",
                             "segments": {"0": {"b1": {'x': 23, 'y': 345.2},
                                                "b2": None,
                                                "from_node_id":"1",
                                                "to_node_id":"2"},
                                      }
                             }},
         "nodes": {"0": {"connected_segments": [{"reaction_id": "0",
                                                 "segment_id": "3"}],
                         "x":3343.8,
                         "y":2167.23,
                         "node_is_primary": True,
                         "label_x":3373.8,
                         "label_y":2177.2,
                         "name": "D-Glucose",
                         "bigg_id": "glc__D_e",
                         "node_type": "metabolite"}
               },
         "text_labels": {},
         "membranes": {},
         "info": {
             "max_map_w": 100,
             "max_map_h": 234.9 }
    }
    conforms_to_spec(o)

    o['nodes']["0"]['x'] = "astring"
    with raises(Exception):
        conforms_to_spec(o)

    o['nodes']["0"]['x'] = None
    with raises(Exception):
        conforms_to_spec(o)

    o['nodes']['x'] = 10
    o['nodes'][3] = o['nodes']['x']
    with raises(Exception):
        conforms_to_spec(o)
        
    del o['nodes']
    with raises(Exception):
        conforms_to_spec(o)
