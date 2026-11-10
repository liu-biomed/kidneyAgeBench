import json, urllib.request, urllib.parse
q = '(ABSTRACT:"aging-related genes" OR ABSTRACT:"ageing-related genes" OR ABSTRACT:"senescence-related genes" OR ABSTRACT:"aging-related differentially expressed genes" OR ABSTRACT:"aging related genes" OR ABSTRACT:"cellular senescence-related" OR ABSTRACT:"ARGs" OR ABSTRACT:"SRGs" OR ABSTRACT:"ARDEGs") AND (TITLE:"kidney" OR TITLE:"renal" OR TITLE:"nephropathy" OR TITLE:"nephritis" OR TITLE:"glomerul*" OR TITLE:"IgA" OR TITLE:"lupus" OR TITLE:"FSGS" OR TITLE:"podocyt*" OR TITLE:"transplant*" OR TITLE:"dialysis" OR TITLE:"CKD")'
out = {}; cur = "*"
for page in range(5):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND SRC:MED", "format": "json", "pageSize": 200, "cursorMark": cur, "resultType": "core"})
    d = json.load(urllib.request.urlopen(u, timeout=90)); rs = d["resultList"]["result"]
    for r in rs: out[r.get("pmid")] = dict(pmid=r.get("pmid"), pmcid=r.get("pmcid"), year=r.get("pubYear"), journal=(r.get("journalInfo") or {}).get("journal", {}).get("title"), title=r.get("title"), oa=r.get("isOpenAccess"), abstract=r.get("abstractText") or "")
    if not rs or d.get("nextCursorMark") in (None, cur): break
    cur = d["nextCursorMark"]
json.dump(list(out.values()), open("/tmp/epmcC.json", "w")); print(len(out))
for v in sorted(out.values(), key=lambda x: x["year"] or ""): print(v["year"], "|", v["pmid"], "|", v["pmcid"], "|", (v["journal"] or "")[:25], "|", v["title"][:110])
