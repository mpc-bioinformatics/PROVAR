##
##    Code ist für das erstellen der tabellen da
##



#!/usr/bin/python3
from configparser import ConfigParser
import psycopg2
from psycopg2 import pool
 

## Connection to the right database
def establish_global_db_connection():
    threaded_postgreSQL_pool = psycopg2.pool.ThreadedConnectionPool(1, 20,user = "postgres",password = "quality",host = "127.0.0.1", port = "5432", database = "provar_backup")
    # Use getconn() method to Get Connection from connection pool
    conn  = threaded_postgreSQL_pool.getconn()
    cur = conn.cursor()
    return conn, cur, threaded_postgreSQL_pool



# ### Connection to the right database
# def establish_global_db_connection():
#     threaded_postgreSQL_pool = psycopg2.pool.ThreadedConnectionPool(1, 20,user = "provar",password = "f4srb2FONOQiG2n",host = "127.0.0.1", port = "5432", database = "provar")
#     # Use getconn() method to Get Connection from connection pool
#     conn  = threaded_postgreSQL_pool.getconn()
#     cur = conn.cursor()
#     return conn, cur, threaded_postgreSQL_pool





def closing_db_connection():
    cur.close()
    #Use this method to release the connection object and send back ti connection pool
    threaded_postgreSQL_pool.putconn(conn)
    #closing database connection.
    # use closeall method to close all the active connection if you want to turn of the application
    if (threaded_postgreSQL_pool):
        threaded_postgreSQL_pool.closeall()

    
def create_tables(conn, cur):
    """ create tables in the PostgreSQL database"""
    commands = (
        """
         CREATE TABLE Anzahl_varianten_pro_protein (
            protein_ID BIGSERIAL PRIMARY KEY , 
            accession VARCHAR(20), 
            anzahl INT,
            len_aa_seq_ INT,
            ratio DECIMAL,
            in_DB INT,
            UNIQUE(accession)
         )
        """,
        """
        CREATE TABLE peptide (
            peptide_ID BIGSERIAL PRIMARY KEY , 
            weight BIGINT, 
            sequence VARCHAR(50) NOT NULL,
            UNIQUE(sequence) 
        );
        """,
        """ CREATE TABLE protein (
                protein_ID BIGSERIAL PRIMARY KEY, 
                accession VARCHAR(20),
                UNIQUE(accession)
                );
        """,
        """
        CREATE TABLE peptide_protein(
                protein_ID BIGSERIAL ,
                peptide_ID BIGSERIAL, 
                organism_ID BIGSERIAL, 
                mutation TEXT,
                start_position INT,
                end_position INT,
                PRIMARY KEY (protein_ID, peptide_ID)
                );
        """,  
        """
        CREATE TABLE organism(
                organism_id BIGSERIAL ,
                organism_name TEXT,
                PRIMARY KEY(organism_name)
                );
        """,
        """
            CREATE TABLE duration_table(
                id BIGSERIAL,
                accession VARCHAR(20),
                number_variants INT,
                duration DECIMAL, 
                start_time DECIMAL,
                end_time DECIMAL
            );
        """
        """
            CREATE TABLE blotting_table(
                time BIGSERIAL,
                count_peptide_protein BIGSERIAL,
                count_peptide BIGSERIAL, 
                number_of_parsed_peptides BIGSERIAL
            );
        """
        )
    
    # create table one by one
    for command in commands:
        cur.execute(command)
    conn.commit()

 
def drop_tables(conn, cur):
    """ create tables in the PostgreSQL database"""
    commands = (
        #"""
        #CREATE TABLE Anzahl_varianten_pro_protein (
        #    protein_ID BIGSERIAL PRIMARY KEY , 
        #    accession VARCHAR(20), 
        #    anzahl INT,
        #    UNIQUE(accession)
        #)
        #""",
        """
        Drop Table peptide;
        """,
        """ Drop Table protein;
        """,
        """
            Drop table peptide_protein;
        """,  
        """
            Drop table organism;
        """,
        """
            Drop table duration_table;
        """,
        """ 
            Drop table blotting_table;
        """
        )
    
    # create table one by one
    for command in commands:
        cur.execute(command)
    conn.commit()    




(conn, cur, threaded_postgreSQL_pool) = establish_global_db_connection()
try:
    drop_tables(conn,cur)
    create_tables(conn,cur)
except:
    create_tables(conn, cur)


closing_db_connection()