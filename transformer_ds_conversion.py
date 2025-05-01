import pandas as pd
import requests
import time

k = 30
INPUT_CSV = "whole_dataset.csv"
OUTPUT_CSV = f"mutation_{k}mer_dataset153.csv"
ENSEMBL_REST = "https://rest.ensembl.org"


def fetch_reference_seq(chrom, start, end, strand=1):
    url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{start}..{end}:{strand}"
    headers = {"Content-Type": "text/plain"}
    print("url ", url)
    response = requests.get(url, headers=headers)
    if not response.ok:
        response.raise_for_status()
    return response.text.strip().upper()


def get_kmer_around_mutation(pos, chrom, ref_nt, alt_nt, k=30, strand=1):
    half_k = k // 2
    start = pos - half_k
    end = pos + half_k - 1
    ref_seq = fetch_reference_seq(chrom, start, end, strand)
    if not ref_seq or len(ref_seq) != k:
        return None, None

    mut_seq = list(ref_seq)
    mut_seq[half_k] = alt_nt.upper()
    return ref_seq, ''.join(mut_seq)


def fetch_protein_sequence(transcript_id):
    url = f"{ENSEMBL_REST}/sequence/id/{transcript_id}?type=protein"
    headers = {"Content-Type": "application/json"}
    r = requests.get(url, headers=headers)
    if r.ok:
        return r.json()["seq"]
    return None


def fetch_vep_info(hgvs):
    url = f"{ENSEMBL_REST}/vep/human/hgvs/{hgvs}"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    r = requests.get(url, headers=headers)
    print(url)
    if r.ok:
        return r.json()
    return None


def process_row(row):
    transcript = row["Transcript"]
    coding = row["c. Change"]
    gene = row["Gene"]
    label = 1 if "deleterious" in row["SIFT Prediction"] else 0

    hgvs = f"{transcript}:{coding}"
    vep_data = fetch_vep_info(hgvs)
    if not vep_data:
        print("fuck")
        return None

    tx_conseq = None
    for cons in vep_data[0].get("transcript_consequences", []):
        if cons.get("transcript_id") == transcript and "sift_score" in cons:
            tx_conseq = cons
            break

    if not tx_conseq:
        return None

    pos = tx_conseq.get("protein_start")
    ref_aa, alt_aa = tx_conseq.get("amino_acids", "/").split("/")
    codons = tx_conseq.get("codons", "///").split("/")
    codon_wt, codon_mut = codons[0].upper(), codons[1].upper()

    prot_seq = fetch_protein_sequence(transcript)
    if not prot_seq or pos > len(prot_seq):
        return None

    start = pos - (k // 2)
    end = pos + (k // 2) - 1
    wt_kmer = prot_seq[start - 1:end]
    if len(wt_kmer) != k:
        return None
    mut_kmer = wt_kmer[:k//2] + alt_aa + wt_kmer[k//2 + 1:]

    chr = vep_data[0].get("seq_region_name")
    genome_pos = vep_data[0].get("start")
    allele_string = vep_data[0].get("allele_string")
    ref_base, alt_base = allele_string.split("/")

    wt_genome_kmer, mut_genome_kmer = get_kmer_around_mutation(genome_pos, chr, ref_base, alt_base, k)

    if not wt_genome_kmer or not mut_genome_kmer:
        return None

    return {
        "Gene": gene,
        "Transcript": transcript,
        "c. Change": coding,
        "Position": pos,
        "WildTypeAA": ref_aa,
        "MutantAA": alt_aa,
        "Codon_WT": codon_wt,
        "Codon_Mut": codon_mut,
        f"WT_{k}mer": wt_kmer,
        f"Mut_{k}mer": mut_kmer,
        f"WT_genome_{k}mer": wt_genome_kmer,
        f"Mut_genome_{k}mer": mut_genome_kmer,
        "Label": label
    }


def main():
    df = pd.read_csv(INPUT_CSV)
    results = []

    for i, row in df.iterrows():
        try:
            data = process_row(row)
            if data:
                results.append(data)
        except Exception as e:
            print(f"Error on row {i}: {e}")

        out_df = pd.DataFrame(results)
        out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()