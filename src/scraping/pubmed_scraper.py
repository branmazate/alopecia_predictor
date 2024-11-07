import os

from dotenv import load_dotenv
from Bio import Entrez, Medline
import pandas as pd
import sqlite3

# Setting up the email address (required by PubMed)
load_dotenv()
email = os.getenv("EMAIL")
Entrez.email = email

# Search term
search_term = "androgenetic alopecia AND diagnosis AND progression"

#Searching articles' IDs and details
def search_pubmed(query, retmax=100):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    resultados = Entrez.read(handle)
    handle.close()
    return resultados['IdList']

def get_details(id_list):
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), retmode="text", rettype='medline')
    data = handle.read()
    handle.close()
    return data

#Parssing the data
def parse_details(data):
    records = []
    for record in Medline.parse(data.splitLines()):
        #Relevant fields extractedd from the record
            pubmed_id= record.get('PMID', ''),
            title = record.get('TI', ''),
            abstract = record.get('AB', ''),
            autors = record.get('AU', ''),
            pubdate = record.get('DP', '')
            
            records.append({
                'pubmed_id': pubmed_id,
                'title': title,
                'abstract': abstract,
                'autors': autors,
                'publication date': pubdate
            })
    return records

#Ejecution
ids = search_pubmed(search_term, retmax=1000)
details = get_details(ids)
regs = parse_details(details)
df = pd.DataFrame(regs)

#Storing the data
connection = sqlite3.connect('/data/raw/pubmed_data.db')
df.to_sql('pubmed_data', connection, if_exists='replace', index=False)
connection.close()