#!/usr/bin/env python
"""
Database tools.
Requires:
pymongo
You need to start a mongo session first with:
mongod --logpath /Volumes/s/G/database_mongo/log --fork --dbpath /Volumes/s/G/database_mongo/
"""
__updated__ = "2017-06-02"

from pymongo import MongoClient

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
    
if __name__ == "__main__":
    import doctest
    doctest.testmod(verbose=True)
