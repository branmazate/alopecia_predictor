import re
import nltk
import spacy
from nltk.corpus import stopwords

#Download the stopwords
nltk.download('stopwords')
stop_words = set(stopwords.words('english'))

#Language model
nlp = spacy.load('en_core_web_sm')

def clean_text(text):
    # Convert to lowercase
    text = text.lower()
    
    # Remove special characters
    text = re.sub(r'[^a-zA-Z\s]', '', text)
    
    #Remove stopwords
    words = text.split()
    words = [word for word in words if word not in stop_words]
    
    #Lemmatization
    doc = nlp(' '.join(words))
    lemmatized_words = [token.lemma_ for token in doc]
    
    #Remove short words
    final_words = [word for word in lemmatized_words if len(word) > 2]
    
    #Join the words back into a string
    clean_text = " ".join(final_words)
    return clean_text

