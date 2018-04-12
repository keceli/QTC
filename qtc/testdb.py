db= [{'_id'     : '[O][O]', 
'stoich'  : 'O2', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:0, 298:0, 'sigma':0.0},
'g09'     : { 
    'b3lyp'   : {'6-311+g(d,p)':{'energy': -150.26605, 'sigma':0}},
    'm062x'   : {'6-311+g(d,p)':{'energy': -150.308589, 'sigma':0}},
    'ccsd' :    {'cc-pvtz':{'energy': -150.1110736, 'sigma':0}},
    'ccsd(t)' : {'cc-pvtz':{'energy': -150.1290363, 'sigma':0,
                          'zmat':'\n 0   3\n o1\n o2                   o1   r2\nVariables:\n R2\t1.2123\n'},
                          'freq': ' '
                           }
            }
},

{'_id'    : '[H][H]', 
'stoich'  : 'H2', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:0, 298:0, 'sigma':0.0},
'g09'     : { 
     'b3lyp'   : {'6-31+g(d,p)':{'energy': -1.17853935, 'sigma':0}},
     'ccsd' :    {'cc-pvtz':{'energy': -1.1723367, 'sigma':0}},
     'ccsd(t)' : {'cc-pvtz':{'energy': -1.1723367, 'sigma':0,
                        'zmat'  : '\n 0   1\n h1\n h2                   h1   r2\nVariables:\n R2\t0.7427\n'}}
            }
},


{'_id'    : 'O', 
'stoich'  : 'H2O', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-238.931, 298:-241.834, 'sigma':0.027},
'g09'     : { 
    'b3lyp'   : {'6-311+g(d,p)':{'energy': -76.434049, 'sigma':0}},
    'ccsd'    : {'cc-pvtz':{'energy': -76.3245455, 'sigma':0}},
    'ccsd(t)' : {'cc-pvtz':{'energy': -76.3322164, 'sigma':0,
                        'zmat'  : '\n 0   1\n o1\n h2                   o1   r2\n h3                   o1   r3       h2   a3\nVariables:\n R2\t0.9594\n R3\t0.9594\n A3\t103.5862\n'}}
            }
},


{'_id'    : 'C', 
'stoich'  : 'CH4', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-66.55, 298:-74.520, 'sigma':0.057},
'g09'     : { 
    'b3lyp'   : {'6-31+g(d,p)':{'energy': -40.52613767, 'sigma':0,
                                  'zmat': '\n 0   1\n c1\n h2                   c1   R1\n h3                   c1   R2       h2   A2\n h4                   c1   R3       h2   A3       h3   D3       0\n h5                   c1   R4       h2   A4       h3   D4       0\nVariables:\n R1\t1.0922\n R2\t1.0922\n A2\t109.4712\n R3\t1.0922\n A3\t109.4712\n D3\t240.0\n R4\t1.0922\n A4\t109.4712\n D4\t120.001\n'}},
    'm062x'   : {'6-311+g(d,p)':{'energy': -40.4967602, 'sigma':0}},
    'ccsd'    : {'cc-pvtz':{'energy': -40.4318181, 'sigma':0}},
    'ccsd(t)' : {'cc-pvtz':{'energy': -40.4380993, 'sigma':0,
                          'zmat': '\n 0   1\n c1\n h2                   c1   R1\n h3                   c1   R2       h2   A2\n h4                   c1   R3       h2   A3       h3   D3       0\n h5                   c1   R4       h2   A4       h3   D4       0\nVariables:\n R1\t1.0889\n R2\t1.0889\n A2\t109.4712\n R3\t1.0889\n A3\t109.4712\n D3\t240.0\n R4\t1.0889\n A4\t109.4712\n D4\t120.0\n'}}
            }
},

{'_id'    : 'O=C=O', 
'stoich'  : 'CO2', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-393.109, 298:-393.475, 'sigma':0.015},
'g09'     : { 
    'b3lyp'   : {'6-311+g(d,p)':{'energy': -188.59039, 'sigma':0,
                               'zmat': '\n 0   1\n o1\n c2                   o1   r2\n o3                   c2   r3       o1   a3\nVariables:\n R2\t1.1694\n R3\t1.1694\n A3\t179.9963\n'}}
            }
},

{'_id'    : 'CC', 
'stoich'  : 'CH3CH3', 
'charge'  :  0,
'mult'    : 1,
'delHf': {0:-68.29, 'sigma':0.14},
'g09'     : { 
    'b3lyp'   : {'6-311+g(d,p)':{'energy': -79.84164, 'sigma':0,
                                   'zmat': '\n 0   1\n c1\n h2                   c1   R1\n h3                   c1   R2       h2   A2\n h4                   c1   R3       h2   A3       h3   D3       0\n h5                   c1   R4       h2   A4       h3   D4       0\nVariables:\n R1\t1.0889\n R2\t1.0889\n A2\t109.4712\n R3\t1.0889\n A3\t109.4712\n D3\t240.0\n R4\t1.0889\n A4\t109.4712\n D4\t120.0\n'}}
            }
},

{'_id'    : 'CO', 
'stoich'  : 'CH3OH', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-189.83, 'sigma':0.16},
'g09'     : { 
    'b3lyp'   : {'6-311+g(d,p)':{'energy': -115.73487, 'sigma':0}}
            }
},

{'_id'    : 'N', 
'stoich'  : 'NH3', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-38.562,298.:-45.554, 'sigma':0.03}
},

{'_id'    : 'N#N', 
'stoich'  : 'N2', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0.:0., 298.:0., 'sigma':0.0}
},

{'_id'    : 'C=C', 
'stoich'  : 'C2H6', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:60.96, 'sigma':0.13}
},

{'_id'    : 'C=O', 
'stoich'  : 'CH2O', 
'charge'  : 0,
'mult'    : 1,
'delHf': {0:-105.349, 298:-109.188, 'sigma':0.099},
'g09'     : { 
    'm062x'   : {'6-311+g(d,p)':{'energy': -114.51152, 'sigma':0,
                                   'zmat':'\n 0   1\n c1\n o2                   c1   R1\n h3                   c1   R2       o2   A2\n h4                   c1   R3       o2   A3       h3   D3       0\nVariables:\n R1\t1.1969\n R2\t1.1047\n A2\t121.7284\n R3\t1.1047\n A3\t121.7284\n D3\t179.9995\n'}}
            }
}]
