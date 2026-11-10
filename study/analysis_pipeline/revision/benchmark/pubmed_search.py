import json, urllib.request, urllib.parse, time, sys, subprocess
def post(url, data):
    for i in range(5):
        r = subprocess.run(['curl','-sS','--max-time','300','--data-binary','@-',url], input=urllib.parse.urlencode(data).encode(), capture_output=True)
        if r.returncode == 0 and r.stdout: return r.stdout.decode()
        time.sleep(3)
    raise SystemExit('curl failed '+url)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
Q = {
 "C": '((aging[tiab] OR ageing[tiab] OR senescence[tiab] OR senescent[tiab] OR "aging-related"[tiab] OR "ageing-related"[tiab] OR "senescence-related"[tiab] OR "senescence-associated"[tiab] OR SASP[tiab]) '
      'AND (kidney[tiab] OR renal[tiab] OR nephropathy[tiab] OR nephritis[tiab] OR glomerul*[tiab] OR podocyte*[tiab] OR "chronic kidney disease"[tiab] OR CKD[tiab] OR DKD[tiab] OR "acute kidney injury"[tiab] OR "kidney transplant*"[tiab] OR "IgA nephropathy"[tiab] OR "lupus nephritis"[tiab] OR "renal cell carcinoma"[tiab] OR ccRCC[tiab] OR KIRC[tiab]) '
      'AND (bioinformatic*[tiab] OR "machine learning"[tiab] OR "hub gene*"[tiab] OR WGCNA[tiab] OR LASSO[tiab] OR "Gene Expression Omnibus"[tiab] OR "GEO database"[tiab] OR TCGA[tiab] OR "prognostic signature"[tiab] OR "risk model"[tiab] OR "diagnostic model"[tiab]) '
      'AND (gene*[tiab] OR signature*[tiab] OR biomarker*[tiab])) AND ("2010/01/01"[dp] : "2026/10/31"[dp]) AND english[la]',
 "A": '((aging[tiab] OR ageing[tiab] OR "donor age"[tiab] OR "chronological age"[tiab]) AND (kidney[tiab] OR renal[tiab]) AND (human[tiab] OR humans[mh]) '
      'AND (transcriptom*[tiab] OR microarray*[tiab] OR "gene expression profil*"[tiab] OR "RNA-seq"[tiab] OR "RNA sequencing"[tiab] OR "expression signature*"[tiab] OR proteom*[tiab]) '
      'AND (age-associated[tiab] OR age-related[tiab] OR "with age"[tiab] OR "aging kidney"[tiab] OR "kidney aging"[tiab] OR "renal aging"[tiab] OR "ageing kidney"[tiab] OR "kidney ageing"[tiab])) AND english[la]'}
out = {}
for k, q in Q.items():
    d = json.loads(post(E + "esearch.fcgi", {"db": "pubmed", "term": q, "retmax": 5000, "retmode": "json"}))
    ids = d["esearchresult"]["idlist"]; print(k, "count", d["esearchresult"]["count"], file=sys.stderr)
    recs = []
    for i in range(0, len(ids), 200):
        x = post(E + "efetch.fcgi", {"db": "pubmed", "id": ",".join(ids[i:i+200]), "rettype": "xml"})
        recs.append(x); time.sleep(0.4)
    out[k] = {"query": q, "count": d["esearchresult"]["count"], "date": time.strftime("%Y-%m-%d"), "ids": ids, "xml": recs}
json.dump(out, open("pubmed_search.json", "w")); print("done")
