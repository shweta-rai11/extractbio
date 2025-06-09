import os
import re
import time
import requests
import pandas as pd
import sqlite3
from Bio import Entrez
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor
import streamlit as st
from datetime import datetime, date
import http.client
from requests.exceptions import ConnectionError

Entrez.email = "swetaraibms@gmail.com"

st.title("Search Keyword and Extract PubMed IDs Linked to Clinical Trial IDs")

def sanitize_keyword(keyword):
    return re.sub(r'[^A-Za-z0-9_]', '_', keyword.strip())

def search_geo_accessions(keyword, retmax=10, retstart=0, retries=3, delay=2):
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="gds", term=keyword, retmax=retmax, retstart=retstart)
            record = Entrez.read(handle)
            handle.close()
            return record.get("IdList", []), int(record.get("Count", 0))
        except RuntimeError as e:
            st.warning(f"NCBI search failed (attempt {attempt+1}/{retries}): {e}")
            time.sleep(delay)
        except Exception as e:
            st.error(f"Unexpected error on search_geo_accessions: {e}")
            break
    return [], 0


def fetch_geo_accession_details(id_list, retries=3):
    accession_numbers = set()
    if id_list:
        for attempt in range(retries):
            try:
                handle = Entrez.efetch(db="gds", id=id_list, rettype="full", retmode="text")
                data = handle.read()
                handle.close()
                pattern = r"GSE\d{1,10}"
                found_accessions = re.findall(pattern, data)
                accession_numbers.update(found_accessions)
                break
            except http.client.IncompleteRead:
                st.warning(f"IncompleteRead on attempt {attempt + 1}/{retries}. Retrying...")
                time.sleep(2 ** attempt)
            except Exception as e:
                st.error(f"❌ Error fetching accession details: {str(e)}")
                break
    return list(accession_numbers)

def fetch_all_geo_accessions(keyword, max_results=100):
    all_accessions = set()
    retstart = 0
    total_count = None
    pbar = st.progress(0)
    while True:
        geo_ids, total_count = search_geo_accessions(keyword, retmax=max_results, retstart=retstart)
        if not geo_ids:
            break
        accessions = fetch_geo_accession_details(geo_ids)
        all_accessions.update(accessions)
        retstart += max_results
        if retstart >= total_count:
            break
        pbar.progress(min(retstart / total_count, 1.0))
    pbar.empty()
    return list(all_accessions)

def fetch_pubmed_ids_from_geo(accession_list):
    base_url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

    def fetch_pubmed_data(accession):
        url = f"{base_url}?acc={accession}&view=full"
        retries = 3
        for attempt in range(retries):
            try:
                response = requests.get(url, timeout=10)
                response.raise_for_status()
                break
            except (ConnectionError, requests.exceptions.RequestException):
                time.sleep(2 ** attempt)
                if attempt == retries - 1:
                    return {'accession': accession, 'Pubmed_ID': ''}
        soup = BeautifulSoup(response.text, 'html.parser')
        pubmed_ids = [re.search(r'pubmed/(\d+)', a['href']).group(1)
                       for a in soup.find_all('a', href=True) if 'pubmed' in a['href'] and re.search(r'pubmed/(\d+)', a['href'])]
        return {'accession': accession, 'Pubmed_ID': ', '.join(pubmed_ids)}

    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()

    with ThreadPoolExecutor(max_workers=10) as executor:
        for i, result in enumerate(executor.map(fetch_pubmed_data, accession_list), 1):
            results.append(result)
            progress_bar.progress(i / len(accession_list))
            status_text.text(f"Processed {i}/{len(accession_list)}: {result['accession']}")

    progress_bar.empty()
    status_text.empty()
    return results

def fetch_pubmed_html(pubmed_id):
    url = f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}/"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        return response.text
    except:
        return None

def get_publication_date(html):
    soup = BeautifulSoup(html, 'html.parser')
    pub_date_div = soup.find('span', class_='cit')
    if pub_date_div:
        text = pub_date_div.get_text()
        match = re.search(r'\d{4} \w{3} \d{1,2}', text) or re.search(r'\d{4}', text)
        if match:
            try:
                return datetime.strptime(match.group(), "%Y %b %d")
            except:
                try:
                    return datetime.strptime(match.group(), "%Y")
                except:
                    return None
    return None

def search_nct_in_abstract(html):
    soup = BeautifulSoup(html, 'html.parser')
    text = soup.get_text(" ", strip=True)
    match = re.search(r'NCT\d{8}', text)
    return match.group(0) if match else None

def process_pubmed_ids(metadata, batch_size=50):
    pubmed_ids = sorted(set(pid.strip() for entry in metadata if entry['Pubmed_ID']
                            for pid in entry['Pubmed_ID'].split(',')))
    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()
    total = len(pubmed_ids)

    for i in range(0, total, batch_size):
        batch = pubmed_ids[i:i+batch_size]
        for j, pid in enumerate(batch):
            html = fetch_pubmed_html(pid)
            nct = search_nct_in_abstract(html) if html else 'NCT Not Found'
            pub_date = get_publication_date(html) if html else None
            results.append({
                'Pubmed_ID': pid,
                'NCT Number': nct if nct else 'NCT Not Found',
                'Publication_Date': pub_date.date() if isinstance(pub_date, datetime) else None
            })
            time.sleep(0.25)
        progress = (i + len(batch)) / total
        progress_bar.progress(min(progress, 1.0))
        status_text.text(f"Processed {i + len(batch)} of {total} PubMed IDs")

    progress_bar.empty()
    status_text.empty()
    return results

def filter_nct(df, start_date, end_date):
    df["Publication_Date"] = pd.to_datetime(df["Publication_Date"], errors='coerce').dt.date
    mask = (df["Publication_Date"].notna()) & \
           (df["Publication_Date"] >= start_date) & \
           (df["Publication_Date"] <= end_date)
    return df.loc[mask & df['NCT Number'].str.match(r'^NCT\d+$', na=False)]

# === STREAMLIT UI ===
keyword = st.text_input("Enter keyword for GEO search:")
start_date = st.date_input("Start publication date", key='start_date')
end_date = st.date_input("End publication date", key='end_date')

if st.button("\U0001F680 Run Search"):
    if not keyword:
        st.error("❗ Please enter a keyword to search.")
    elif not start_date or not end_date:
        st.error("❗ Please select both start and end publication dates.")
    else:
        with st.spinner("Fetching GEO accessions..."):
            accessions = fetch_all_geo_accessions(keyword)
        geo_df = pd.DataFrame(accessions, columns=["GEO Accession Number"])
        st.dataframe(geo_df)

        with st.spinner("Fetching PubMed metadata..."):
            metadata = fetch_pubmed_ids_from_geo(accessions)
        metadata_df = pd.DataFrame(metadata)
        st.dataframe(metadata_df)

        with st.spinner("Extracting NCT Numbers & Publication Dates..."):
            nct_df = pd.DataFrame(process_pubmed_ids(metadata))
        st.dataframe(nct_df)

        with st.spinner("Filtering based on publication date..."):
            filtered_nct_df = filter_nct(nct_df, start_date, end_date)
        st.dataframe(filtered_nct_df)

        # Save all data
        db_name = "combined_data.db"
        conn = sqlite3.connect(db_name)
        geo_df.to_sql("geo_data", conn, if_exists='replace', index=False)
        metadata_df.to_sql("pubmed_data", conn, if_exists='replace', index=False)
        nct_df.to_sql("nct_data", conn, if_exists='replace', index=False)
        conn.commit()
        conn.close()

        st.success("✅ All data saved to one combined SQLite database.")

        with open(db_name, "rb") as f:
            st.download_button(
                label="\U0001F4E5 Download Combined DB",
                data=f,
                file_name=db_name,
                mime="application/x-sqlite3"
            )
