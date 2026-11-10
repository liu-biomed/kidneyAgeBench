import json, urllib.request, urllib.parse, sys
qs = {
 "orig": '(TITLE:"kidney" OR TITLE:"renal" OR TITLE:"nephron") AND (TITLE:"aging" OR TITLE:"ageing" OR TITLE:"age" OR TITLE:"aged" OR TITLE:"older" OR TITLE:"elderly") AND (TITLE:"transcript*" OR TITLE:"gene expression" OR TITLE:"microarray" OR TITLE:"RNA" OR TITLE:"molecular signature*" OR TITLE:"genes")',
 "hub": '(ABSTRACT:"kidney aging" OR ABSTRACT:"renal aging" OR ABSTRACT:"kidney ageing" OR ABSTRACT:"renal ageing" OR ABSTRACT:"aging kidney" OR ABSTRACT:"aged kidney") AND (ABSTRACT:"GSE*" OR ABSTRACT:"Gene Expression Omnibus" OR ABSTRACT:"differentially expressed genes" OR ABSTRACT:"hub genes" OR ABSTRACT:"WGCNA")',
 "sc": '(ABSTRACT:"kidney aging" OR ABSTRACT:"renal aging" OR ABSTRACT:"aging kidney" OR ABSTRACT:"kidney ageing" OR ABSTRACT:"aging human kidney") AND (ABSTRACT:"single-cell" OR ABSTRACT:"single-nucleus" OR ABSTRACT:"scRNA-seq" OR ABSTRACT:"snRNA-seq" OR ABSTRACT:"spatial transcriptom*")',
 "clock": '(ABSTRACT:"kidney" OR ABSTRACT:"renal") AND (ABSTRACT:"transcriptomic age" OR ABSTRACT:"transcriptomic clock" OR ABSTRACT:"RNA age" OR ABSTRACT:"aging clock" OR ABSTRACT:"organ age" OR ABSTRACT:"tissue-specific aging signature*") AND (ABSTRACT:"GTEx" OR ABSTRACT:"transcriptom*" OR ABSTRACT:"proteom*")',
 "donor": '(ABSTRACT:"living donor*" OR ABSTRACT:"zero-hour" OR ABSTRACT:"0-hour" OR ABSTRACT:"pre-implantation" OR ABSTRACT:"procurement biops*") AND (ABSTRACT:"donor age") AND (ABSTRACT:"gene expression" OR ABSTRACT:"transcriptom*" OR ABSTRACT:"microarray")',
}
allh = {}
for k, q in qs.items():
    cur = "*"; n = 0
    for page in range(6):
        u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND SRC:MED", "format": "json", "pageSize": 200, "cursorMark": cur, "resultType": "core"})
        d = json.load(urllib.request.urlopen(u, timeout=90)); rs = d["resultList"]["result"]
        for r in rs:
            key = r.get("pmid") or r.get("id"); n += 1
            if key not in allh: allh[key] = dict(cat=k, pmid=r.get("pmid"), pmcid=r.get("pmcid"), doi=r.get("doi"), year=r.get("pubYear"), journal=(r.get("journalInfo") or {}).get("journal", {}).get("title"), title=r.get("title"), abstract=(r.get("abstractText") or "")[:2500], cited=r.get("citedByCount"))
            else: allh[key]["cat"] += "," + k
        if not rs or d.get("nextCursorMark") in (None, cur): break
        cur = d["nextCursorMark"]
    print(k, n, file=sys.stderr)
json.dump(list(allh.values()), open("/tmp/epmc2.json", "w")); print("unique", len(allh))
