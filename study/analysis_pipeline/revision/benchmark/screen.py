import json, re, csv, html
d = json.load(open("/Users/zheranliu/kidney_aging_data/lit/pubmed_search.json")); sym = set(l.strip() for l in open("/tmp/hgnc_pc.txt"))
bad_sym = {"AGE","AKI","CKD","DKD","GEO","TCGA","ROC","AUC","PCR","LASSO","SVM","RF","DN","LN","SLE","MCD","FSGS","GSEA","PPI","KEGG","IGA","NK","CI","HR","OS","SD","CS","SRG","ARG","MR","II","III","IV","PD","TME","ICI","RCC","MAP","CAT","SET","REST","PC","CA","GC","LPS","STAR","ACE","RAGE","SASP","ESRD","IGAN","HE","IHC","MIR","LIGHT","IMPACT","ROS","AR","CD","T","B","GO","DE","DEG","DEGS","NR","MS","CR","EPO","PH","ALB","TIMP"}
rows = []; seen = set()
for k in ["C", "A"]:
    for x in d[k]["xml"]:
        for art in re.findall(r"<PubmedArticle>.*?</PubmedArticle>", x, re.S):
            g = lambda p: (re.search(p, art, re.S) or [None, ""])[1]
            pmid = g(r"<PMID[^>]*>(\d+)</PMID>")
            if pmid in seen: continue
            seen.add(pmid)
            title = html.unescape(re.sub(r"<[^>]+>", "", g(r"<ArticleTitle>(.*?)</ArticleTitle>")))
            ab = html.unescape(re.sub(r"<[^>]+>", " ", " ".join(re.findall(r"<AbstractText[^>]*>(.*?)</AbstractText>", art, re.S))))
            pt = re.findall(r"<PublicationType[^>]*>(.*?)</PublicationType>", art)
            yr = g(r"<PubDate>.*?<Year>(\d{4})</Year>"); jn = g(r"<ISOAbbreviation>(.*?)</ISOAbbreviation>"); pmc = g(r'<ArticleId IdType="pmc">(PMC\d+)</ArticleId>')
            toks = re.findall(r"\b[A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]+)?\b", ab); genes = [t for t in dict.fromkeys(toks) if t in sym and t not in bad_sym]
            t = (title + " " + ab).lower()
            review = any(p in ("Review", "Systematic Review", "Meta-Analysis", "Editorial", "Comment") for p in pt)
            aging = bool(re.search(r"aging|ageing|senescen|sasp", t)); kidney = bool(re.search(r"kidney|renal|nephr|glomerul|podocyt|kirc|ccrcc|dialysis", t))
            model = bool(re.search(r"hub gene|signature|lasso|machine learning|wgcna|risk model|prognostic model|diagnostic model|biomarker|feature gene|key gene|characteristic gene", t))
            animal_only = bool(re.search(r"\b(mice|mouse|rat|rats|murine)\b", t)) and not re.search(r"\bhuman|patient", t)
            mirna = bool(re.search(r"lncrna|mirna|circrna|long non-coding|microrna", title.lower()))
            other = bool(re.search(r"autophagy|anoikis|cuproptosis|ferroptosis|pyroptosis|necroptosis|disulfidptosis|m6a|glycolysis|lactylation", title.lower())) and not re.search(r"aging|ageing|senescen", title.lower())
            rows.append(dict(search=k, pmid=pmid, pmcid=pmc, year=yr, journal=jn, title=title, review=review, aging=aging, kidney=kidney, model=model, animal_only=animal_only, ncRNA=mirna, other_death=other, n_genes_abs=len(genes), genes_abs=";".join(genes), abstract=ab))
out = "/Users/zheranliu/kidney_aging_data/lit/screen_all.csv"
with open(out, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
inc = [r for r in rows if r["search"] == "C" and not r["review"] and r["aging"] and r["kidney"] and r["model"] and not r["animal_only"] and not r["ncRNA"] and not r["other_death"]]
print("records", len(rows), "| C-search", sum(r["search"] == "C" for r in rows), "| auto-eligible C", len(inc))
for r in inc: print(r["pmid"], r["pmcid"] or "-", r["year"], r["n_genes_abs"], "|", r["title"][:95], "|", r["genes_abs"][:90])
