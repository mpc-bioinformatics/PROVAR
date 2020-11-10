##
##   OCde ist für das erstellen der Fasta-Datei da
##

import time
import sys
from datetime import datetime
import numpy as np
from psycopg2 import pool
import psycopg2
from multiprocessing import Process, Event, Queue, Pool
from multiprocessing.queues import Empty as EmptyQueueError


### Connection to the right database
def establish_global_db_connection():
    threaded_postgreSQL_pool = psycopg2.pool.ThreadedConnectionPool(1, 20,user = "postgres",password = "quality",host = "127.0.0.1", port = "5432", database = "test_database")
    # Use getconn() method to Get Connection from connection pool
    conn  = threaded_postgreSQL_pool.getconn()
    cur = conn.cursor()
    return conn, cur, threaded_postgreSQL_pool



def closing_db_connection(conn, cur, threaded_postgreSQL_pool):
    cur.close()
    #Use this method to release the connection object and send back ti connection poolf4srb2FONOQiG2n
    threaded_postgreSQL_pool.putconn(conn)
    #closing database connection.
    # use closeall method to close all the active connection if you want to turn of the application
    if (threaded_postgreSQL_pool):
        threaded_postgreSQL_pool.closeall()



def get_info(peptide_id,conn,cur):
    cmd = """select 
            (protein.accession, peptide_protein.mutation,peptide_protein.start_position, peptide_protein.end_position, peptide.weight, peptide.sequence)
            from peptide
            join peptide_protein on peptide_protein.peptide_id = peptide.peptide_id
            join organism on organism.organism_id = peptide_protein.organism_id
            join protein on peptide_protein.protein_id = protein.protein_id
            where peptide_protein.peptide_id = """+str(peptide_id)+"""
            order by peptide.peptide_id"""
    cur.execute(cmd)
    protein_ids = cur.fetchall()
    return protein_ids


def organisms(organism_id, conn ,cur):
    cmd = """SELECT peptide_id 
            from peptide_protein 
            WHERE organism_id = """+str(organism_id)+""" 
            group by peptide_id;
            """
    cur.execute(cmd)
    protein_ids = cur.fetchall()
    return protein_ids


def get_organism(SQL_organism_ID, conn, cur):
    cmd = "Select organism_name from organism WHERE organism_id = "+str(SQL_organism_ID)+";"
    cur.execute(cmd)
    organism = cur.fetchone()[0]
    return organism
def how_many_organisms_available(conn, cur):
    cmd = "Select Count(*) from organism;"
    cur.execute(cmd)
    available_organisms = cur.fetchone()[0]
    return available_organisms  

def worker(fasta_queue, stop_flag: Event, queue: Queue):
    (conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()
    while not (stop_flag.is_set() and queue.empty()):                    ###!!!!!!!!
        peptide_ids = queue.get()       
        while len(peptide_ids) != 0:
            peptide_id = peptide_ids.pop(0)[0]
            output = get_info(peptide_id, conn,cur)[0][0]
            output = output[1:-1] 
            output = output.split(",")
            fasta_info = [output[0] ,output[1][1:-1] ,output[2] ,output[3],output[4] ,output[5] ]      #Accession, mutation, start_pos, end_pos, MW, sequence <-- Reihenfolge, wie print in fasta soll
            fasta_queue.put(fasta_info)
        else:
            print("finished process", datetime.now())

    closing_db_connection(conn, cur, threaded_postgreSQL_pool)



def writer(save_file ,stop_flag: Event, queue: Queue):
    while not (stop_flag.is_set() and queue.empty()): 
        write_this = fasta_queue.get()
        save_file.write(">"+write_this[0] +":"+write_this[1]+"_positions:"+str(write_this[2])+"-"+str(write_this[3])+"_MW:"+str(write_this[4])+"Da \n"+write_this[5]+ "\n")
    else:
        print("finished writing", datetime.now())


def create_indices(conn, cur):
    cmd = "CREATE INDEX peptide_protein_protein_id_idx ON peptide_protein (protein_ID);"
    cmd2 = "CREATE INDEX peptide_protein_peptide_id_idx ON peptide_protein (peptide_ID);"
    cmd3 = "CREATE INDEX peptide_protein_organism_id_idx ON peptide_protein (organism_ID);"
    try:
        cur.execute(cmd)
        conn.commit()
    except:
        pass
    try:
        cur.execute(cmd2)
        conn.commit()
    except:
        pass    
    try:
        cur.execute(cmd3)
    except:
        pass
    print("Indices sind erstellt")

def chunks(lst,n):
    return [lst[i::n] for i in range(n)]
 
print("Start:", datetime.now() )
(conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()   #Verbindung
#create_indices(conn, cur)
available_organisms = how_many_organisms_available(conn, cur)               
organism_id = 41                        # nur für diesen Code --> Mensch an position 41 in DB
organism_name = get_organism(organism_id, conn, cur)    #für dieses mal irrelevent, sonst wichtig für fasta-file bezeichnung
save_file=open('/home/baranme/Documents/Master/output/'+'HomoSapiens_30M_output(2).fasta', 'w')     
peptide_ids = organisms(organism_id, conn, cur)      #Alle wichtigen ids, die abgearbeitet/abgefragt werden sollen
print(len(peptide_ids), "are available")

peptide_chunks = chunks(peptide_ids, 11)       #Ähnliche Verteilung für die Prozesse erstellen, die eingelesen werden
print("splitted into", len(peptide_chunks), "even parts")

closing_db_connection(conn, cur, threaded_postgreSQL_pool) #wird ab hier nicht mehr benötigt

fasta_queue = Queue()
worker_queue = Queue()

fasta_stop_flag = Event()
worker_stop_flag= Event()


while len(peptide_chunks) != 0:
    worker_queue.put(peptide_chunks.pop(0))                 #einlesen in Queue
print("was given to the processes", datetime.now())

procs = []
for _ in range(11):   #ABLESEPROZESSE
    new_proc = Process(target=worker, args=(fasta_queue, worker_stop_flag, worker_queue,)) # The last ',' in args is necessary!
    new_proc.start()
    procs.append(new_proc)

second_procs = []
for idx in range(1):   #1 Weil wenn zwei, dann halten die sich gegenseitig auf!            #EINLESEPROZESS
    second_new_proc = Process(target=writer, args=(save_file, fasta_stop_flag ,fasta_queue,))
    second_new_proc.start()
    second_procs.append(second_new_proc)

worker_stop_flag.set()
for proc in procs:
    proc.join() 
    
fasta_stop_flag.set()
for second_proc in second_procs:
    second_proc.join()


save_file.close()
