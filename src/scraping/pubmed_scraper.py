import os

from dotenv import load_dotenv
from Bio import Entrez, Medline
from io import StringIO
import pandas as pd
import sqlite3

# Setting up the email address (required by PubMed)
load_dotenv()
email = os.getenv("EMAIL")
Entrez.email = email

# Search term
search_term = "androgenetic alopecia AND diagnostic AND progression"

# Searching articles' IDs and details
def search_pubmed(query, retmax=100):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    resultados = Entrez.read(handle)
    handle.close()
    print(f'Found {resultados["Count"]} articles')
    return resultados['IdList']

def get_details(id_list):
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list),
                        retmode="text", rettype='medline')
    data = handle.read()
    handle.close()
    return data

# Parssing the data
def parse_details(data):
    handle = StringIO(data)
    
    records = []
    for record in Medline.parse(handle):
        #Relevant fields extractedd from the record
            pubmed_id= record.get('PMID', '')
            title = record.get('TI', '')
            abstract = record.get('AB', '')
            authors = record.get('AU', '')
            pubdate = record.get('DP', '')

            records.append({
                'pubmed_id': pubmed_id,
                'title': title,
                'abstract': abstract,
                'authors': authors,
                'publication date': pubdate
            })
    handle.close()
    return records

# Ejecution
ids = search_pubmed(search_term, retmax=1000)
details = get_details(ids)
regs = parse_details(details)
df = pd.DataFrame(regs)

#Turning authors in a string
df['authors'] = df['authors'].apply(lambda x: ', '.join(x) if isinstance(x, list) else x)

# Storing the data in a sqlite database
connection = sqlite3.connect('data/raw/pubmed_data.db') 

df.to_sql('pubmed_data', connection, if_exists='replace', index=False)

connection.close()
