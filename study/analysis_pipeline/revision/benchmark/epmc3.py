import json, urllib.request, urllib.parse
qs = ['(TITLE:"renal aging" OR TITLE:"kidney aging" OR TITLE:"renal ageing" OR TITLE:"kidney ageing" OR TITLE:"aging kidney" OR TITLE:"aged kidney") AND (ABSTRACT:"bioinformatic*" OR ABSTRACT:"GEO" OR ABSTRACT:"differentially expressed" OR ABSTRACT:"hub" OR ABSTRACT:"machine learning" OR ABSTRACT:"signature")',
      '(TITLE:"human kidney" OR TITLE:"human kidneys") AND (TITLE:"aging" OR TITLE:"ageing" OR TITLE:"age") AND (ABSTRACT:"single-cell" OR ABSTRACT:"single-nucleus" OR ABSTRACT:"transcriptom*" OR ABSTRACT:"proteom*")',
      '(ABSTRACT:"tissue-specific" OR ABSTRACT:"across tissues" OR ABSTRACT:"multi-tissue" OR ABSTRACT:"pan-tissue") AND ABSTRACT:"GTEx" AND (TITLE:"aging" OR TITLE:"ageing" OR TITLE:"age") AND (ABSTRACT:"kidney")',
      '(TITLE:"kidney" OR TITLE:"renal") AND (TITLE:"proteom*") AND (TITLE:"aging" OR TITLE:"ageing" OR TITLE:"age") AND ABSTRACT:"human"']
seen = {}
for q in qs:
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND SRC:MED", "format": "json", "pageSize": 200, "resultType": "lite"})
    for r in json.load(urllib.request.urlopen(u, timeout=90))["resultList"]["result"]: seen.setdefault(r.get("pmid"), (q[:12], r.get("pubYear"), r.get("pmid"), r.get("pmcid"), (r.get("journalTitle") or "")[:26], r.get("title")[:125]))
for v in sorted(seen.values(), key=lambda x: (x[0], x[1])): print(" | ".join(str(z) for z in v))
