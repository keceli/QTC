#!/usr/bin/env python
"""
Database tools.
Requires:
pymongo
You need to start a mongo session first with:
mongod --logpath /Volumes/s/G/database_mongo/log --fork --dbpath /Volumes/s/G/database_mongo/
"""
__updated__ = "2017-06-23"
try:
    from pymongo import MongoClient
except:
    pass
import json
import logging
def start_mongo(dbpath,logpath,mongoexe='mongod'):
    """
    Starts a mongo session.
    
    """
    import subprocess
    subprocess.call([mongoexe,'--logpath', logpath, '--fork', '--dbpath', dbpath])
    return


def get_collection(c='mycollection'):
    from pymongo import MongoClient
    client = MongoClient()
    return client[c]


def get_database(collection,db='mydatabase'):
    return collection[db]


def insert_entry(entry,db,efilter={}):
    return db.replace_one(efilter,entry,upsert=True)


def get_smiles_from_json(jsonfile):
    """
    Opens and reads a .json file created by RMG. 
    Builds a list of strings contains SMILES.
    Needs to convert from unicode to str.
    """
    with open(jsonfile) as f:
        j = json.load(f)
    nitem = len(j)
    smileslist = ['']*nitem
    i = 0
    for d in j :
        smileslist[i] = d['SMILES'].encode('ascii','ignore')
        i += 1
    return smileslist

def gen_dict_extract(key, var):
    """
    https://stackoverflow.com/questions/9807634/find-all-occurrences-of-a-key-in-nested-python-dictionaries-and-lists
    """
    if hasattr(var,'iteritems'):
        for k, v in var.iteritems():
            if k == key:
                yield v
            if isinstance(v, dict):
                for result in gen_dict_extract(key, v):
                    yield result
            elif isinstance(v, list):
                for d in v:
                    for result in gen_dict_extract(key, d):
                        yield result


def visualise_dict(d,lvl=0):
    """
    https://stackoverflow.com/questions/15023333/simple-tool-library-to-visualize-huge-python-dict
    """
    # go through the dictionary alphabetically 
    for k in sorted(d):
        form = '{:<45} {:<15} {:<10}'
        # print the table header if we're at the beginning
        if lvl == 0 and k == sorted(d)[0]:
            logging.info(form.format('KEY','LEVEL','TYPE'))
            logging.info('-'*79)

        indent = '  '*lvl # indent the table to visualise hierarchy
        t = str(type(d[k]))

        # print details of each entry
        logging.info(form.format(indent+str(k),lvl,t))

        # if the entry is a dictionary
        if type(d[k])==dict:
            # visualise THAT dictionary with +1 indent
            visualise_dict(d[k],lvl+1)
    
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
