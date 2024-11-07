import os

from dotenv import load_dotenv
from Bio import Entrez

# Setting up the email address (required by PubMed)
load_dotenv()
email = os.getenv("EMAIL")
Entrez.email = email

# Search term
search_term = "androgenetic alopecia AND diagnosis AND progression"

#Searching articles' IDs 
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