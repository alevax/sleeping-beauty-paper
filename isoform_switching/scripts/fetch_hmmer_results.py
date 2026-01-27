#!/usr/bin/env python3
"""
Fetch HMMER results from API and convert to webserver email format for IsoformSwitchAnalyzeR.

Expected format (tab-delimited, 16 columns):
seq id, alignment start, alignment end, envelope start, envelope end,
hmm acc, hmm name, type, hmm start, hmm end, hmm length, bit score,
E-value, significance, clan
"""

import json
import requests
import sys
import time

def fetch_result(job_id):
    """Fetch a single job result."""
    url = f"https://www.ebi.ac.uk/Tools/hmmer/api/v1/result/{job_id}"
    response = requests.get(url, headers={"Accept": "application/json"})
    return response.json()

def main(jobs_file, output_file):
    # Load job list
    with open(jobs_file) as f:
        jobs = json.load(f)

    print(f"Fetching {len(jobs)} results...")

    # Open output file and write header
    with open(output_file, 'w') as out:
        # Header row - this is what the webserver sends via email
        out.write("seq id\talignment start\talignment end\tenvelope start\tenvelope end\thmm acc\thmm name\ttype\thmm start\thmm end\thmm length\tbit score\tE-value\tsignificance\tclan\n")

        for i, job in enumerate(jobs):
            job_id = job['id']
            query_name = job['query_name']

            try:
                result = fetch_result(job_id)

                if result.get('status') != 'SUCCESS':
                    print(f"  {query_name}: status {result.get('status')}")
                    continue

                hits = result.get('result', {}).get('hits', [])

                if not hits:
                    print(f"  {query_name}: no hits")
                    continue

                for hit in hits:
                    if not hit.get('is_included', False):
                        continue

                    metadata = hit.get('metadata', {})
                    hmm_name = metadata.get('identifier', 'NA')
                    hmm_acc = metadata.get('accession', 'NA')
                    hmm_type = metadata.get('type', 'Domain')
                    hmm_length = metadata.get('model_length', 0)
                    clan = metadata.get('clan', '')
                    if clan is None:
                        clan = ''

                    domains = hit.get('domains', [])

                    for domain in domains:
                        if not domain.get('is_included', False):
                            continue

                        # Extract domain info
                        dom_score = domain.get('bitscore', 0)
                        dom_evalue = domain.get('ievalue', 1)

                        align = domain.get('alignment_display', {})
                        hmm_from = align.get('hmmfrom', 0)
                        hmm_to = align.get('hmmto', 0)
                        ali_from = align.get('sqfrom', 0)
                        ali_to = align.get('sqto', 0)
                        env_from = domain.get('ienv', 0)
                        env_to = domain.get('jenv', 0)

                        # Significance (1 = significant)
                        significance = 1

                        # Format: seq_id, ali_start, ali_end, env_start, env_end, hmm_acc, hmm_name, type, hmm_start, hmm_end, hmm_length, bit_score, e-value, significance, clan
                        line = f"{query_name}\t{ali_from}\t{ali_to}\t{env_from}\t{env_to}\t{hmm_acc}\t{hmm_name}\t{hmm_type}\t{hmm_from}\t{hmm_to}\t{hmm_length}\t{dom_score:.1f}\t{dom_evalue:.2e}\t{significance}\t{clan}\n"
                        out.write(line)

                print(f"  [{i+1}/{len(jobs)}] {query_name}: {len(hits)} hits")

            except Exception as e:
                print(f"  {query_name}: ERROR - {e}")

            # Small delay to be nice to the API
            time.sleep(0.1)

    print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    jobs_file = sys.argv[1] if len(sys.argv) > 1 else "results/isoformswitch/cis_overlap/hmmer_jobs.json"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "results/isoformswitch/cis_overlap/pfam_results.txt"
    main(jobs_file, output_file)
