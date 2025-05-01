import pandas as pd
import re
import requests
from time import sleep

INPUT_FILE = "datasets\GeneDatasheet.xlsx"
BATCH_SIZE = 200


def extract_variant_parts(name_str):
    try:
        refseq = re.search(r"(NM_\d+(?:\.\d+)?)", name_str).group(1)
        gene = re.search(r"\(([^)]+)\)", name_str).group(1)
        coding = re.search(r"c\.\-?\d+[ACGT]>[ACGT]", name_str).group(0)
        return refseq, gene, coding
    
    except:
        return None, None, None


def strip_refseq_version(nm_id):
    return nm_id.split(".")[0]


def get_ensembl_transcript(nm_id):
    base_id = strip_refseq_version(nm_id)
    url = f"https://rest.ensembl.org/xrefs/name/homo_sapiens/{base_id}?content-type=application/json"
    r = requests.get(url)
    if r.ok:
        for item in r.json():
            if item.get("id", "").startswith("ENST"):
                return item["id"]
    return None


def get_all_transcripts_by_gene(gene):
    url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?expand=1;content-type=application/json"
    r = requests.get(url)
    if r.ok:
        return [t['id'] for t in r.json().get('Transcript', []) if t['id'].startswith("ENST")]
    return []


def batch_get_sift_from_vep(hgvs_list):
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json" }
    r = requests.post(vep_url, headers=headers, json={"hgvs_notations": hgvs_list})
    if not r.ok:
        print(f"VEP batch request failed: {r.status_code}")
        return {}

    results = {}
    try:
        for variant in r.json():
            hgvs = variant["input"]
            for t in variant.get("transcript_consequences", []):
                if "sift_score" in t:
                    results[hgvs] = {
                        "sift_score": t["sift_score"],
                        "sift_prediction": t["sift_prediction"]
                    }
                    break
    except Exception as e:
        print("Batch parse error:", e)
    return results


def main():
    df = pd.read_excel(INPUT_FILE)
    df = df[df["Molecular consequence"].str.contains("missense", case=False, na=False)]

    all_variants = []
    metadata_map = {}
    seen_keys = set()

    for _, row in df.iterrows():
        name = row["Name"]
        refseq, gene, coding = extract_variant_parts(name)
        if not coding:
            continue

        key = (refseq, coding)
        if key in seen_keys:
            continue
        seen_keys.add(key)

        direct_tx = get_ensembl_transcript(refseq)
        transcripts = [direct_tx] if direct_tx else get_all_transcripts_by_gene(gene)

        for tx in transcripts:
            hgvs = f"{tx}:{coding}"
            all_variants.append(hgvs)
            metadata_map[hgvs] = {
                "Gene": gene,
                "RefSeq": refseq,
                "Transcript": tx,
                "c. Change": coding
            }
            break  # only try the first transcript for batching
    
    results = []
    for i in range(0, len(all_variants), BATCH_SIZE):
        batch = all_variants[i:i + BATCH_SIZE]
        print(f"⏱ Processing batch {i//BATCH_SIZE + 1} with {len(batch)} variants...")
        batch_result = batch_get_sift_from_vep(batch)

        for hgvs in batch:
            meta = metadata_map[hgvs]
            if hgvs in batch_result:
                r = batch_result[hgvs]
                print(f"✅ {hgvs}: {r['sift_prediction']}, {r['sift_score']}")
                results.append({
                    **meta,
                    "SIFT Prediction": r["sift_prediction"],
                    "SIFT Score": r["sift_score"]
                })
            else:
                print(f"no SIFT result for {hgvs}")
        # sleep(1)

    out_df = pd.DataFrame(results)
    out_df.to_csv("sift_predictions_vep_batch.csv", index=False)
    print("saved to sift_predictions_vep_batch.csv")


if __name__ == "__main__":
    main()
