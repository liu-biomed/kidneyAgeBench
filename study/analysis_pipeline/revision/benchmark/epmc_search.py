import json, urllib.request, urllib.parse, sys
qs = [
 '(TITLE:"kidney" OR TITLE:"renal") AND (TITLE:"aging" OR TITLE:"ageing" OR TITLE:"senescence" OR TITLE:"age-related" OR TITLE:"elderly") AND (ABSTRACT:"bioinformatic*" OR ABSTRACT:"GEO" OR ABSTRACT:"Gene Expression Omnibus" OR ABSTRACT:"hub gene*" OR ABSTRACT:"WGCNA") AND (ABSTRACT:"human")',
 '(TITLE:"kidney" OR TITLE:"renal") AND (TITLE:"aging" OR TITLE:"ageing") AND (ABSTRACT:"signature" OR ABSTRACT:"biomarker*" OR ABSTRACT:"hub gene*" OR ABSTRACT:"machine learning")',
 '(TITLE:"kidney" OR TITLE:"renal") AND (TITLE:"aging" OR TITLE:"ageing") AND (ABSTRACT:"transcriptom*" OR ABSTRACT:"microarray" OR ABSTRACT:"RNA-seq" OR ABSTRACT:"RNA sequencing") AND ABSTRACT:"human"',
 '(TITLE:"kidney" OR TITLE:"renal") AND (TITLE:"senescence-related" OR TITLE:"aging-related" OR TITLE:"ageing-related" OR TITLE:"senescence related") AND (TITLE:"gene*" OR ABSTRACT:"signature")',
 '"kidney aging" AND ("transcriptomic clock" OR "aging clock" OR "biological age") AND (transcriptom* OR "gene expression")',
]
seen = {}
for q in qs:
    cur = "*"
    for page in range(4):
        u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND SRC:MED", "format": "json", "pageSize": 100, "cursorMark": cur, "resultType": "lite"})
        d = json.load(urllib.request.urlopen(u, timeout=60)); rs = d["resultList"]["result"]
        for r in rs: seen.setdefault(r.get("pmid") or r.get("id"), r)
        if not rs or d.get("nextCursorMark") in (None, cur): break
        cur = d["nextCursorMark"]
    print("query hits so far", len(seen), file=sys.stderr)
out = [dict(pmid=r.get("pmid"), pmcid=r.get("pmcid"), doi=r.get("doi"), year=r.get("pubYear"), journal=r.get("journalTitle"), title=r.get("title"), cited=r.get("citedByCount")) for r in seen.values()]
json.dump(out, open("epmc_hits.json", "w"), indent=0); print(len(out))
