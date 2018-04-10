#!/usr/local/bin/python3.4
import sqlite3 as lite
from sys import argv

conn = lite.connect('aPrimerCustomDB.db')
cur = conn.cursor()

def get_posts():
    exStr='SELECT * FROM primers WHERE length between '+str(argv[1])+' and '+str(argv[2])
    cur.execute(exStr)
    return(cur.fetchall())

tabs=get_posts()
for i in tabs:
    print(i[2]+"\t"+i[3])
