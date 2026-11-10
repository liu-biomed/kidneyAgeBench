import json, urllib.request, urllib.parse, re
q = '(TITLE:"aging" OR TITLE:"ageing" OR TITLE:"senescence" OR TITLE:"senescent" OR TITLE:"aging-related" OR TITLE:"senescence-related") AND (TITLE:"kidney" OR TITLE:"renal" OR TITLE:"nephropathy" OR TITLE:"nephritis" OR TITLE:"glomerul*" OR TITLE:"podocyte*" OR TITLE:"CKD" OR TITLE:"DKD" OR TITLE:"AKI" OR TITLE:"transplant*") AND (ABSTRACT:"bioinformatic*" OR ABSTRACT:"machine learning" OR ABSTRACT:"hub gene*" OR ABSTRACT:"WGCNA" OR ABSTRACT:"LASSO" OR ABSTRACT:"Gene Expression Omnibus" OR ABSTRACT:"GEO database")'
out = {}; cur = "*"
for page in range(5):
    u = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q + " AND SRC:MED", "format": "json", "pageSize": 200, "cursorMark": cur, "resultType": "core"})
    d = json.load(urllib.request.urlopen(u, timeout=90)); rs = d["resultList"]["result"]
    for r in rs: out[r.get("pmid")] = dict(pmid=r.get("pmid"), pmcid=r.get("pmcid"), year=r.get("pubYear"), journal=(r.get("journalInfo") or {}).get("journal", {}).get("title"), title=r.get("title"), abstract=re.sub(r"<[^>]+>", "", r.get("abstractText") or ""))
    if not rs or d.get("nextCursorMark") in (None, cur): break
    cur = d["nextCursorMark"]
json.dump(list(out.values()), open("/tmp/epmcC2.json", "w")); print(len(out))
have = {"37760894","37396184","37957539","39602449","40671964","40593989","40319112","40290492","42152390","35646056","35720278","34207247","42102994","42227369","42792788","42464365"}
for v in sorted(out.values(), key=lambda x: x["year"] or ""):
    if v["pmid"] in have: continue
    print(v["year"], "|", v["pmid"], "|", v["pmcid"], "|", (v["journal"] or "")[:22], "|", v["title"][:115])
