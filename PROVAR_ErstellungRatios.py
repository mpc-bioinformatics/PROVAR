##
## Diese Code ist für das erstellen der Ratios da
##


import matplotlib
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import write_dot
from multiprocessing import Process
from threading import Thread
import os
import math
import time
import cProfile
import regex
import re
from datetime import datetime
from itertools import chain, combinations
import itertools
import math
from configparser import ConfigParser
import psycopg2
from psycopg2 import pool
import networkx as nx
import pstats
import io
import cProfile
import pstats
import io
from pstats import SortKey
from multiprocessing import Process, Event, Queue, Pool
from multiprocessing.queues import Empty as EmptyQueueError
from concurrent.futures import ProcessPoolExecutor
import concurrent


def filtering_info(data):
    #aa_seq filtering
    data1 = data.split(";")[-1]
    data1 = data1.replace("//", "")
    data1 = data1.replace("     ","")
    aa_seq = data1.replace(" ", "")
    del(data1)
    #filtering variants
    data2 = data.replace("Missing", "  -> *")   #um die deletierten mutationen zu filtern, muss diese processiert werden
    filtering_only_variant_text = regex.findall(r"(FT   VARIANT         \d*.*.*)", data2)

    variant_char = regex.findall(r"FT   VARIANT         \d*FT.{19}/note=.{6}(.)", data2, regex.MULTILINE)
    variant_coordinaten = regex.findall(r"FT   VARIANT         (\d*)FT.{19}/note=.{6}[\w|*]", data2, regex.MULTILINE)
    variants = [] 
    i = 0
    while i < len(variant_char):
        try:
            eintrag = "("+str(variant_char[i])+"|"+str(variant_coordinaten[i])+")"
            variants.append(str(eintrag))
            i = i + 1
        except:
            break    

    #accession filtering
    accession = regex.findall(r"ID   \w*.*.\nAC   (\w*)", data)
    organism = ''.join(regex.findall(';OS.{3}([^\s]+\s+[^\s]+)', data))
    if '.OC' in organism:
        organism = organism[:-3] 
    if "'" in organism:
        organism = organism.replace("'", "")  
    return accession[0], variants, variant_coordinaten, variant_char, aa_seq, organism


### Connection to the right database
def establish_global_db_connection():
    threaded_postgreSQL_pool = psycopg2.pool.ThreadedConnectionPool(1, 20,user = "provar",password = "f4srb2FONOQiG2n",host = "127.0.0.1", port = "5432", database = "provar")
    # Use getconn() method to Get Connection from connection pool
    conn  = threaded_postgreSQL_pool.getconn()
    cur = conn.cursor()
    return conn, cur, threaded_postgreSQL_pool



def closing_db_connection(conn, cur, threaded_postgreSQL_pool):
    cur.close()
    threaded_postgreSQL_pool.putconn(conn)
    if (threaded_postgreSQL_pool):
        threaded_postgreSQL_pool.closeall()



def count_fkt(accession, conn, cur):
    cmd = """select 
            Count(*)
            from peptide
            join peptide_protein on peptide_protein.peptide_id = peptide.peptide_id
            join protein on peptide_protein.protein_id = protein.protein_id
            where accession = '"""+str(accession)+"""';"""
    cur.execute(cmd)
    wie_viele_eintraege_in_diesem_protein=cur.fetchone()[0] 
    return wie_viele_eintraege_in_diesem_protein



def import_variantanzahl_into_db(accession, anzahl_varianten, aa_seq, ratio ,in_DB, conn, cur):
        cmd = "INSERT INTO anzahl_varianten_pro_protein (protein_id,accession, anzahl, len_aa_seq_, ratio, in_DB) VALUES (DEFAULT,'" +str(accession)+"', " +str(anzahl_varianten)+", "+str(len(aa_seq))+", "+str(ratio)+", "+str(in_DB)+");"
        commands = (
                   [cmd]
                    
                   )

        for command in commands:
            cur.execute(str(command))


protein_acc = ""
protein_seq = ""
(conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()
q = 0
with open('/mnt/data/provar/uniprot_sprot.dat', 'r') as f:
    for line in f:
        if (line.startswith('ID')):
            #altes protein abarbeiten
            if (protein_acc != ""):
                data = protein_acc + "\n" + protein_seq
                accession, variants, variant_coordinaten, variant_char, aa_seq, organism = filtering_info(data)
                anzahl_eintraege = count_fkt(accession, conn, cur)
                ratio = len(variants)/len(aa_seq)
                import_variantanzahl_into_db(accession, len(variants), aa_seq, ratio, anzahl_eintraege, conn, cur)

            protein_acc = line.rstrip('\n')
            protein_seq = ""
            q = q+1
            if q % 1000 ==0:
                conn.commit()
               

        else:
            protein_seq += line.rstrip('\n')

closing_db_connection(conn, cur, threaded_postgreSQL_pool)
