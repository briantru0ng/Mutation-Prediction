import pandas as pd
import re
import requests

# INPUT_FILE = "sift-data.xlsx"
INPUT_FILE = "datasets\GeneDatasheet.xlsx"


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


def get_sift_from_vep(hgvs_string):
    vep_url = "https://rest.ensembl.org/vep/human/hgvs"
    headers = { "Content-Type" : "application/json", "Accept" : "application/json" }
    data = { "hgvs_notations": [hgvs_string] }

    r = requests.post(vep_url, headers=headers, json=data)
    if not r.ok:
        return None

    try:
        result = r.json()[0]
        for t in result.get("transcript_consequences", []):
            if "sift_score" in t:
                return {
                    "sift_score": t["sift_score"],
                    "sift_prediction": t["sift_prediction"]
                }
    except:
        return None

    return None


def main():
    df = pd.read_excel(INPUT_FILE)
    # missense variants only
    df = df[df["Molecular consequence"].str.contains("missense", case=False, na=False)]

    results = []

    for _, row in df.iterrows():
        name = row["Name"]
        refseq, gene, coding = extract_variant_parts(name)
        if not coding:
            continue

        transcripts = []
        direct = get_ensembl_transcript(refseq)
        if direct:
            transcripts.append(direct)
        else:
            transcripts = get_all_transcripts_by_gene(gene)

        found = False
        for tx in transcripts:
            hgvs = f"{tx}:{coding}"
            sift = get_sift_from_vep(hgvs)
            if sift:
                print(f"✅ Found SIFT for {hgvs} ({sift['sift_prediction']}, {sift['sift_score']})")
                results.append({
                    "Gene": gene,
                    "RefSeq": refseq,
                    "Transcript": tx,
                    "c. Change": coding,
                    "SIFT Prediction": sift["sift_prediction"],
                    "SIFT Score": sift["sift_score"]
                })
                found = True
                break

        if not found:
            print(f"no SIFT for {gene}:{coding}")

    out_df = pd.DataFrame(results)
    out_df.to_csv("genedataset_vep.csv", index=False)
    print("done")


if __name__ == "__main__":
    main()
