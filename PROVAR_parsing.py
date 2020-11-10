##
## Code ist für das verdauen, sowie das einschreiben da
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
import concurrent.futures
import psutil

print(psutil.cpu_percent())

weights = {'A': 71.03711, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259, 'F': 147.06841,
           'G': 57.02146, 'H': 137.05891, 'I': 113.08406, 'K': 128.09496, 'L': 113.08406,
           'M': 131.04049, 'N': 114.04293, 'P': 97.05276, 'Q': 128.05858, 'R': 156.10111,
           'S': 87.03203, 'T': 101.04768, 'V': 99.06841, 'W': 186.07931, 'Y': 163.06333 }


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
    threaded_postgreSQL_pool = psycopg2.pool.ThreadedConnectionPool(1, 20, user = "provar",password ="f4srb2FONOQiG2n" , host = "max-decoy-db", port = "5432", database = "provar")
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


###Importing the Proteinnames to the DB (Tablename protein)    
def import_protein_into_db(accession, conn, cur):
   # ###########################################################################
   # ..................execute statements.......................................
   # ###########################################################################     
    
    cmd = "INSERT INTO protein (protein_id, accession) VALUES (DEFAULT,'" +str(accession)+"') RETURNING protein_id;"
    

    
    cur.execute(cmd)
    SQL_protein_ID = cur.fetchone()[0]
    conn.commit()    

    # ###########################################################################
    # ..................end execute statements...................................
    # ###########################################################################  
    return SQL_protein_ID





def getting_inserted_count(conn,cur):
    cmd = "select Count(*) from peptide_protein;"
    cmd2 = "select Count(*) from peptide"

    cur.execute(cmd)
    count_peptide_protein = cur.fetchone()[0]  
    cur.execute(cmd2)
    count_peptide = cur.fetchone()[0]   
    return count_peptide_protein, count_peptide

### Function, um die Varianten auf Peptidebene zu erhalten
def peptide_variants(aa_seq, variants):
    #tryptic digestion
    peptides_after_digestion = regex.split(r'(?<=[RK])(?=[^P])', aa_seq)
    
    ###
    ###Falls es annotierte Varianten existieren, dann tu diese Funktion
    ###
    try:
        variant_coordinates = re.findall("(\d+)", str(variants)) #Indexe, an der die Mutationen geschehen

        j= 0
        summationsliste= []   ###summationsliste = Länge der jeweiligen peptide nach dem Verdau
        while j != len(peptides_after_digestion):
            x = len(peptides_after_digestion[j])
            summationsliste.append(x)
            j = j + 1
                
        zwischensumme_list = []
        j = 1
        while j <= len(summationsliste):   #Finden der zu wechselnden AS durch hochzählen (Zahl kommt von der Datenbank)
            zwischensumme = sum(summationsliste[0:j])
            zwischensumme_list.append(zwischensumme)
            j = j + 1

        #### Darstellen bzw. Erhalten der Varianten auf PEPTIDEBENE
        peptide_variants = [[]for x in range(len(peptides_after_digestion))]
        i = 0
        j = 0
        k = 0
        while k <= len(peptides_after_digestion):
            try:
                while int(variant_coordinates[i]) <= int(zwischensumme_list[j]):
                    peptide_variants[j].append(str(variants[i]))
                    i = i + 1
                j = j + 1
                k = k + 1
                
            except:
                    break
        return peptide_variants, zwischensumme_list    
    except:
        pass    


def import_combinations_info(SQL_peptide_ID, SQL_protein_ID, SQL_organism_ID,variantliste, start_position, end_position, conn, cur):   
    cmd =  "INSERT INTO peptide_protein(protein_ID, peptide_ID,organism_ID, mutation, start_position, end_position) VALUES ("+str(SQL_protein_ID)+" ," +str(SQL_peptide_ID) + ", "+str(SQL_organism_ID) + ", '"+ str(variantliste)+"', "+ str(start_position)+", "+ str(end_position)+") ON CONFLICT DO NOTHING;"
   
    cur.execute(cmd)

###Importing the generated Mutations to the DB (Tablename peptide)    
def import_mutated_peptide_into_db(aa_sequence, molecular_weight, conn, cur):
    cmd = "INSERT INTO peptide (peptide_id,weight, sequence) VALUES (DEFAULT," +str(molecular_weight)+",'"+str(aa_sequence)+"') ON CONFLICT DO NOTHING RETURNING peptide_id ;"
    cur.execute(cmd)
    result = cur.fetchone() 
    if result is None:
        SQL_peptide_ID = check_if_peptide_exists(aa_sequence, conn, cur)
        return SQL_peptide_ID    
    else:
        return result[0]    

def blotting_information(time,count_peptide_protein,count_peptide,anzahl_geparste_processe,conn,cur)   :
    cmd = "INSERT INTO blotting_table(time, count_peptide_protein, count_peptide, number_of_parsed_peptides) VALUES (" +str(time)+", "+str(count_peptide_protein)+", "+str(count_peptide)+", "+str(anzahl_geparste_processe)+");"
    cur.execute(cmd)


def check_if_organism_exists(organism, conn, cur):
    cmd = "SELECT organism_id from organism WHERE organism_name = '"+str(organism)+"';"
    cur.execute(cmd)
    organism_status = cur.fetchone()
    if organism_status:
        return organism_status[0]
    else:
        return None 


def get_accession(SQL_protein_ID, conn, cur):
    cmd = "Select accession from protein WHERE protein_id = "+str(SQL_protein_ID)+";"
    cur.execute(cmd)
    accession = cur.fetchone()[0]
    return accession




def get_organism(SQL_organism_ID, conn, cur):
    cmd = "Select organism_name from organism WHERE organism_id = "+str(SQL_organism_ID)+";"
    cur.execute(cmd)
    organism = cur.fetchone()[0]
    return organism

def import_organism_DB(organism, conn, cur):
    cmd = "INSERT INTO organism (organism_id, organism_name) VALUES (DEFAULT, '" +str(organism)+"') ON CONFLICT DO NOTHING RETURNING organism_id ;"
   
    cur.execute(cmd)
    SQL_organism_ID = cur.fetchone()[0]
    #print("line 157")
    return SQL_organism_ID    


# add the variants-nodes
def adding_variant_edges(seqGraph, variant_coordinaten, variant_char, variants):
    i = 0
    for i in range(0, len(variant_coordinaten)):  
        a = int(variant_coordinaten[i])-1
        b = int(variant_coordinaten[i])
        c = variant_char[i]
        d = variants[i]
        seqGraph.add_edge(a,b, sequence = c, variants = d)
        #print(a,b,c)
        i = i + 1
    #print(i)    


def calculating_molecular_weight(sequence, start_position):
    
    end_position = start_position + len(sequence) - 1

    sequence = sequence.replace("*", "")      
    #sequence = ''.join(re.split(r'(?<=[RK])(?=[^P])', sequence)[0])
    weight_list = [18.010565]
    for acid in sequence:
        weight = sum(weights[a] for a in acid)
        weight_list.append(weight)
    molecular_weight = sum(weight_list)
    molecular_weight = int(molecular_weight * 1000000) 
    return sequence, molecular_weight, end_position

def check_if_peptide_exists(sequence, conn, cur):
    cmd = "SELECT peptide_id from peptide WHERE sequence = '"+str(sequence)+"';"
    cur.execute(cmd)
    # Wenn das peptide nicht existiert, gibt fetchone() None zurück wenn du darauf [0] aufrufst, dann kommt der Fehler von gerade
    peptide_status = cur.fetchone()
    if peptide_status:
        return peptide_status[0]
    else:
        return None    


def check_if_protein_exists(accession, conn, cur):
    cmd = "SELECT EXISTS(SELECT 1 from protein where accession = '"+str(accession)+"');"
   
    cur.execute(cmd)
    protein_status = cur.fetchone()[0]
    #cur.execute("COMMIT")    
    #print("Line 196")
    return protein_status        

def filter_and_append_to_finalEdges(sequence, molecular_weight, variantliste, SQL_protein_ID,SQL_organism_ID,start_position, end_position, finalEdges):
    coor_var = int(re.findall("\d+", variantliste)[0])

    #Hier wird sicher gestellt, dass nur die Sequenzen gespeichert werden, die auch Varianten beinhalten
    if len(variantliste) != 0 and 6 <= len(sequence) <= 50 and start_position < coor_var < end_position  and sequence not in finalEdges[0]:    
        finalEdges[0][sequence] = len(finalEdges[0])
        finalEdges[1].append(variantliste)
        finalEdges[2].append(SQL_protein_ID)
        finalEdges[3].append(SQL_organism_ID)
        finalEdges[4].append(molecular_weight)
        finalEdges[5].append(start_position)
        finalEdges[6].append(end_position)

        #print(sequence, variantliste, SQL_protein_ID, start_position, end_position)
    return finalEdges

class aminoacid:
    def __init__(self, sequence, variants):
        aminoacid.sequence = sequence
        aminoacid.variants = variants




def get_the_peptide_to_be_parsed(zwischensumme_list, variants_on_peptidelevel, aa_seq, SQL_protein_ID, SQL_organism_ID, variant_coordinaten, variant_char, variants, accession   ):
    process_these = [] 
    b = 0
    while b < len(zwischensumme_list):
        info = [] 
        try:
            if len(variants_on_peptidelevel[b]) > 0:
                #info.append(len(variants_on_peptidelevel[b]))   #Anzahl der Varianten im Peptid
                info.append(zwischensumme_list[b])
                info.append(aa_seq)
                info.append(SQL_organism_ID)
                info.append(SQL_protein_ID)
                info.append(variant_coordinaten)
                info.append(variant_char)
                info.append(variants)
                info.append(len(variants_on_peptidelevel[b]))
                info.append(accession)
                process_these.append(info)
        except:
            pass
        b = b + 1
    return process_these

def call_up_variable(finalEdges, duration_item):
    (conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()
    accession = duration_item[0]
    duration = duration_item[1]
    start = duration_item[2]
    end = duration_item[3]  
    anzahl_varianten_am_peptid = duration_item[4] 
    import_duration_to_DB(accession, duration, anzahl_varianten_am_peptid, start, end, conn, cur)
    #save_file=open('/home/baranme/Documents/Master/output/' + str(organism)+'_output.fasta', 'a')#
    if len(finalEdges[0]) != 0:
        for sequence, i in finalEdges[0].items():   
            #sequence = finalEdges[0][i]
            variantliste = ''.join(finalEdges[1][i])
            SQL_protein_ID = int(finalEdges[2][i])
            SQL_organism_ID = int(finalEdges[3][i])
            molecular_weight = int(finalEdges[4][i])
            start_position =  int(finalEdges[5][i])
            end_position =  int(finalEdges[6][i])
            import_to_DB(str(sequence), variantliste, SQL_protein_ID, SQL_organism_ID, molecular_weight, start_position, end_position, conn, cur)   
            #print(str(sequence), variantliste, SQL_protein_ID, SQL_organism_ID, molecular_weight, start_position, end_position)
            #save_file.write(">"+accession+":"+variantliste+"_positions:"+str(start_position)+"-"+str(+end_position)+"_MW:"+str(molecular_weight)+"Da \n"+sequence+ "\n")
            #print(">"+accession+":"+variantliste+"_positions:"+str(start_position)+"-"+str(+end_position)+"_MW:"+str(molecular_weight)+"Da \n"+sequence+ "\n")
            
            
    conn.commit()
    closing_db_connection(conn, cur, threaded_postgreSQL_pool)

    #save_file.close()

def import_to_DB(sequence, variantliste, SQL_protein_ID, SQL_organism_ID, molecular_weight,start_position, end_position, conn, cur) :
    global peptide_counter
    peptide_status = check_if_peptide_exists(sequence, conn, cur)  #look if peptide exists
    #print("peptide_status",peptide_status)
    peptide_counter = peptide_counter + 1   
    #print(sequence, variantliste, SQL_protein_ID, SQL_organism_ID, molecular_weight,start_position, end_position)
    if peptide_status is not None:
        #Wenn peptidsequenz schon da ist, dann nimm diese Peptid_Id und füge die combinationsinfo in die DB ein
        SQL_peptide_ID = peptide_status
        import_combinations_info(SQL_peptide_ID, SQL_protein_ID,SQL_organism_ID,variantliste, start_position, end_position, conn, cur) 
    else:    
        #Wenn eintrag nicht da, füge Peptidsequenz in die DB und nimm diese Peptid_Id und importiere die combinationsinformation 
        #print("neuer peptid", sequence)
        
        SQL_peptide_ID = import_mutated_peptide_into_db(sequence, molecular_weight, conn, cur)  
        import_combinations_info(SQL_peptide_ID, SQL_protein_ID,SQL_organism_ID,variantliste, start_position, end_position, conn, cur)


def processing(sequence, variantliste, SQL_protein_ID,SQL_organism_ID,finalEdges, start_position, u, d, aa_seq):
    #print(sequence, variantliste,u,d, "all")
    new_aa_seq = aa_seq[0:(start_position)] + sequence + aa_seq[(start_position+49):-1]
    
    variantliste = re.findall("\(.\|\d*\)", variantliste)
    variants_peptide1, leere_variable = peptide_variants(new_aa_seq, variantliste)
    #print(variants_peptide1)
    del(leere_variable)
    sequence = str(re.split(r'(?<=[RK])(?=[^P])', sequence)[0])
    try:
        # um nur die Varianten von dem ersten peptid zu erhalten
        variants_peptide = ''.join([entry for entry in variants_peptide1 if len(entry) != 0][0]) 
        if sequence[-1] == "R" or sequence[-1] == "K"  :
            start_position = start_position + 1
            sequence1, molecular_weight, end_position = calculating_molecular_weight(sequence, start_position)
            finalEdges = filter_and_append_to_finalEdges(sequence1, molecular_weight, variants_peptide, SQL_protein_ID,SQL_organism_ID, start_position,end_position,finalEdges) 
            #conn.commit()
            return finalEdges
    except:
            pass


def iterateAllPathsEdgesUtil(graph, u, d, visited, path, pathEdges, finalEdges, variantEdges, nodeSequence, SQL_protein_ID, SQL_organism_ID, aa_seq):
    start_position = u
    #start_time = float(round(time.time(),6))        
    start = time.perf_counter()
    end_all_info= intern_iterateAllPathsEdgesUtil(graph, u, d, visited, path, pathEdges, finalEdges, variantEdges, nodeSequence, SQL_protein_ID,SQL_organism_ID,aa_seq, start_position)  
    
    end = time.perf_counter()
    duration = end - start
    # duration = float(round(time.time(),6)) - start_time
    # duration = round(duration, 6)
    return end_all_info, duration, start, end


def intern_iterateAllPathsEdgesUtil(graph, u, d, visited, path, pathEdges, finalEdges, variantEdges,nodeSequence, SQL_protein_ID,SQL_organism_ID, aa_seq, start_position): 
    """ helper tool to iterate all paths and get the visited edges """
    # Mark the current node as visited and store in path 
    visited[u] = True
    path.append(u)
   
    # If current vertex is same as destination, then work on the result
    # for now, only the edges string is recorded... but these could also already be send to the DB, if teh array grows too big... or have parallel tasks querying the array
    if u == d:
        sequence = "".join(nodeSequence)
        variantliste = "".join(variantEdges)
        finalEdges = processing(sequence, variantliste, SQL_protein_ID,SQL_organism_ID, finalEdges, start_position, u, d, aa_seq)
        # If current vertex is not destination 
        # Recur for all the vertices adjacent to this vertex
    else:
        for a, b, dat in graph.edges(u, data=True):
            if visited[b] == False:
                try:   
                    next_aa = aminoacid(dat['sequence'], dat['variants']) 
                    pathEdges.append(next_aa)                               #Es liegt nicht bei allen Knoten ein Variant vor!
                except:
                    next_aa = aminoacid(dat['sequence'], "") 
                    pathEdges.append(next_aa)  
                nodeSequence.append(''.join(str(pathEdges[0].sequence)))
                variantEdges.append(''.join(str(pathEdges[0].variants)))
                
                intern_iterateAllPathsEdgesUtil(graph, b, d, visited, path, pathEdges, finalEdges, variantEdges, nodeSequence, SQL_protein_ID,SQL_organism_ID,aa_seq, start_position)
                      
    # Remove current vertex from path[] and mark it as unvisited 
    path.pop()
    if (len(pathEdges) > 0):
        pathEdges.pop()
        variantEdges.pop()
        nodeSequence.pop()

        
    visited[u] = False
    return finalEdges



def iterateAllPathsEdges(list_of_chunks):
    
    
    variant_coordinaten = list_of_chunks[4]
    variant_char = list_of_chunks[5]
    variants = list_of_chunks[6] 
    anzahl_der_varianten = list_of_chunks[7]
  
    
    aa_seq = list_of_chunks[1]
    SQL_organism_ID = int(list_of_chunks[2])
    SQL_protein_ID = int(list_of_chunks[3])
    u = int(list_of_chunks[0])
    d = u + 50
    if d > len(aa_seq):
        d = len(aa_seq)

    seqGraph = generate_network(aa_seq, variant_coordinaten, variant_char, variants)
    
    path = []
    pathEdges = []
    finalEdges = [dict(),[],[],[],[],[],[]] 
    variantEdges = [] 
    nodeSequence = [] 
    visited = [False]*len(seqGraph.nodes)
    start_position = int(u)
    end_all_info, duration, start, end = iterateAllPathsEdgesUtil(seqGraph, u, d, visited, path, pathEdges, finalEdges, variantEdges, nodeSequence, SQL_protein_ID, SQL_organism_ID, aa_seq)
    

    return end_all_info, duration, start, end


    
def chunks(list_elements, n):
   ''' Generate n long list_elements '''
   for i in range(0, len(list_elements), n):
        yield list_elements[i:i+n]





def generate_network(aa_seq, variant_coordinaten, variant_char, variants):
    seqGraph = nx.MultiDiGraph()   
    pos = 0
    seqGraph.add_node(pos)
    for c in aa_seq:
        lastPos = pos
        pos += 1
        seqGraph.add_node(pos)
        seqGraph.add_edge(lastPos, pos, sequence=c)

    adding_variant_edges(seqGraph, variant_coordinaten, variant_char, variants) 

    return seqGraph

def chunks_by_sum1(seq, n):
    mean = sum(map(lambda x:x[0],seq))/n
    seq = sorted(seq, key=lambda  x:(int(x[0]), x[1] ))
    while seq:
        chunk = [seq.pop(-1)]  #[[20, 'a']]
        
        partial_sum = 0
        while seq and partial_sum < mean:
            partial_sum = sum(map(lambda x:x[0], chunk))  #20
            m = len(seq)      #11
            i = 0
            j = 0
            min_diff = abs(partial_sum + seq[i][0] - mean)    #1.0
            last_diff = min_diff
            while i < m and last_diff > min_diff:
                diff = abs(partial_sum + x - mean)
                if diff < min_diff:
                    min_diff = diff
                    j = i
                last_diff = min_diff
                i += 1
            partial_sum += seq[j][0]
            chunk.append(seq.pop(j))
        yield chunk

def import_duration_to_DB(accession, duration, variantenlaenge, start, end, conn,cur):
    cmd = "INSERT INTO duration_table (id, accession, number_variants, duration, start_time, end_time) VALUES (DEFAULT,'" +str(accession)+"',"+str(variantenlaenge)+", "+str(duration)+", "+str(start)+", "+str(end)+") ON CONFLICT DO NOTHING ;"
    cur.execute(cmd)



def parser(SQL_queue,time_worker, stop_flag: Event, queue: Queue):
    # Don't stop until stop flag is set AND queue is empty
    while not (stop_flag.is_set() and queue.empty()):                    ###!!!!!!!!
        try:
            # queue.get(False) returns immediately and throws queue.Empty if queue is empty 
            item = queue.get(False)   
            anzahl_varianten_am_peptid = item[-2]
            accession = item[-1] 

            end_all_info, duration, start, end = iterateAllPathsEdges(item)
            
            info_time_worker = [accession, duration, start, end, anzahl_varianten_am_peptid] 

            #adding to SQL_queue
            end_all_info = end_all_info[0:7] 
            info_for_end= (end_all_info, info_time_worker)
            SQL_queue.put(info_for_end)
            #time_worker.put(info_time_worker)

        except EmptyQueueError:
            # Sleep for 200 ms if the queue is empty
            time.sleep(0.2)
        except KeyboardInterrupt:
            pass




def all_into_DB(stop_flag: Event, queue: Queue, time_worker: Queue):
    while not stop_flag.is_set() or not queue.empty():                  
        try:
            # queue.get(False) returns immediately and throws queue.Empty if queue is empty 
            item = queue.get(timeout=2)  
            sequences_info=item[0]  
            duration_item = item[1] 
            #duration_item = time_worker.get(timeout=2)
            #adding to SQL_queue
            call_up_variable(sequences_info, duration_item)

            

            
        except EmptyQueueError:
            # Sleep for 200 ms if the queue is empty
            pass
        except KeyboardInterrupt:
            pass        

#start profile

protein_acc = ""
protein_seq = ""
(conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()
peptide_counter = 0
worker_queue = Queue()
SQL_queue = Queue()
time_worker = Queue()

parser_stop_flag = Event()
second_consumer_stop_flag = Event()
time_stop_flag = Event

end_all_info = [] 
q = 0      # ist die protein_ID
with open('/mnt/data/provar/uniprot_sprot.dat', 'r') as f:    
    for line in f:
        if (line.startswith('ID')):
            #altes protein abarbeiten
            if (protein_acc != ""):
                data = protein_acc + "\n" + protein_seq
                accession, variants, variant_coordinaten, variant_char, aa_seq, organism = filtering_info(data)
                organism_status = check_if_organism_exists(organism, conn, cur)
                if organism_status is not None:
                    SQL_organism_ID = organism_status 
                else:
                    SQL_organism_ID = import_organism_DB(organism, conn, cur)
                SQL_protein_status = check_if_protein_exists(accession, conn, cur)
                
                if SQL_protein_status == False:
                    
                    SQL_protein_ID = import_protein_into_db(accession, conn, cur)   
                    conn.commit()
                    variants_on_peptidelevel, zwischensumme_list = peptide_variants(aa_seq, variants)
                    zwischensumme_list.insert(0,0)     #Hier die anfangskoordinaten für die Peptide
                    if len(variant_coordinaten) == 0:
                        print(accession, datetime.now() , q)
                    if len(variant_coordinaten) != 0:
                        #Erstellen des Graphen, sowie das dranhängen der Knoten und Kanten
                        process_these_peptides = get_the_peptide_to_be_parsed(zwischensumme_list, variants_on_peptidelevel, aa_seq, SQL_protein_ID, SQL_organism_ID, variant_coordinaten, variant_char, variants, accession)

                        i = 0
                        while i<len(process_these_peptides):
                            worker_queue.put(process_these_peptides[i])
                            i = i+1
                        
                        print(accession, datetime.now() , q, "len var",len(variants), worker_queue.qsize(), "are in the queue")  

                        #process_these_peptides = process_these_peptides[19:23] 
            protein_acc = line.rstrip('\n')
            protein_seq = ""
            q = q+1
            if q % 2 ==0:
                conn.commit()
               

        else:
            protein_seq += line.rstrip('\n')


procs = []

#letztes protein
data = protein_acc + "\n" + protein_seq
accession, variants, variant_coordinaten, variant_char, aa_seq, organism = filtering_info(data)
organism_status = check_if_organism_exists(organism, conn, cur)
if organism_status is not None:
    SQL_organism_ID = organism_status 
else:
    SQL_organism_ID = import_organism_DB(organism, conn, cur)
SQL_protein_status = check_if_protein_exists(accession, conn, cur)
if SQL_protein_status == False:
    
    SQL_protein_ID = import_protein_into_db(accession, conn, cur)   
    conn.commit()
    variants_on_peptidelevel, zwischensumme_list = peptide_variants(aa_seq, variants)
    zwischensumme_list.insert(0,0)     #Hier die anfangskoordinaten für die Peptide
    if len(variant_coordinaten) == 0:
        print(accession, datetime.now() , q)
    if len(variant_coordinaten) != 0:
        #Erstellen des Graphen, sowie das dranhängen der Knoten und Kanten
        process_these_peptides = get_the_peptide_to_be_parsed(zwischensumme_list, variants_on_peptidelevel, aa_seq, SQL_protein_ID, SQL_organism_ID, variant_coordinaten, variant_char, variants, accession)

        i = 0
        while i<len(process_these_peptides):
            worker_queue.put(process_these_peptides[i])
            i = i+1
        
        print(accession, datetime.now() , q, "len var", len(variants), worker_queue.qsize(), "are in the queue")  



def anzahl_der_processe_die_schon_durch_sind(conn, cur):
    cmd = "select Count(*) from duration_table;"

    cur.execute(cmd)
    anzahl_geparste_processe = cur.fetchone()[0]
    return anzahl_geparste_processe  


procs = []
for _ in range(83):
    # Create new proc
    new_proc = Process(target=parser, args=(SQL_queue, time_worker, parser_stop_flag, worker_queue,)) # The last ',' in args is necessary!
    # Start new proc
    new_proc.start()
    # Add proc to the proc array
    procs.append(new_proc)

second_procs = []
for idx in range(7):
    second_new_proc = Process(target=all_into_DB, args=(second_consumer_stop_flag,SQL_queue,time_worker, ))
    second_new_proc.start()
    second_procs.append(second_new_proc)

condition_for_while = 0

counter_time = 0
lenght_of_blot = 900   #Ab welchen abständen soll ein Count gezählt werden
#while worker_queue.qsize() != 0 or SQL_queue.qsize() != 0:
while psutil.cpu_percent() >= 0.1 : 
    count_peptide_protein, count_peptide = getting_inserted_count(conn, cur)
    anzahl_geparste_processe = anzahl_der_processe_die_schon_durch_sind(conn,cur)
    blotting_information(counter_time, count_peptide_protein, count_peptide, anzahl_geparste_processe, conn, cur)
    counter_time = counter_time + lenght_of_blot
    conn.commit()
    time.sleep(lenght_of_blot)

parser_stop_flag.set()

# Wait for each consumer to stop
for proc in procs:
    proc.join() 


# Set the stop flag, so the consumers stop when the queue is empty
second_consumer_stop_flag.set()

# Wait for each consumer to stop
for second_proc in second_procs:
    second_proc.join()

condition_for_while = 1

# lenght_of_schelife = False


conn.commit()
closing_db_connection(conn, cur, threaded_postgreSQL_pool)

print("THE END",datetime.now())


##### making visual graph with gephi
#write_dot(seqGraph, "test.dot")
#print("run ./gephi-0.9.2/bin/gephi and choose your dot file ")

