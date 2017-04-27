#!/usr/bin/env python

"""
Quantum chemistry driver.
"""

def get_args():
    """
    Returns args object that contains command line options.
    """
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=
    """
    April 11, 2017
    Murat Keceli
    
    
    """)
#     parser.add_argument('deltaH', type=float, help='Delta H0 [kcal/mol], float')
#     parser.add_argument('tag',    type=str,   help='Comment tag, string')
#     parser.add_argument('name',   type=str,   help='Species formula, string',nargs='?')
    parser.add_argument('-r','--runserial',action='store_true',nargs=0,
                        help='Run serial')
#     parser.add_argument('-m','--messpf',type=argparse.FileType('r'),nargs=1,
#                         default='/tcghome/ygeorgi/fock/crossrate/bin/partition_function',
#                         help='Path for mess partition function executable')
#     parser.add_argument('-t','--thermp',type=argparse.FileType('r'),nargs=1,
#                         default='/tcghome/sjk/gen/aux_me/therm/thermp.exe',
#                         help='Path for thermp executable')
#     parser.add_argument('-p','--pac99',type=argparse.FileType('r'),nargs=1,
#                         default='/tcghome/sjk/gen/aux_me/therm/pac99.x',
#                         help='Path for pac99 executable')
    return parser.parse_args()