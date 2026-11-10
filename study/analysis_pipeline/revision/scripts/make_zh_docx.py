import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

def set_cell_background(cell, fill_hex):
    tcPr = cell._element.get_or_add_tcPr()
    shd = OxmlElement('w:shd')
    shd.set(qn('w:val'), 'clear')
    shd.set(qn('w:color'), 'auto')
    shd.set(qn('w:fill'), fill_hex)
    tcPr.append(shd)

def set_cell_margins(cell, top=100, bottom=100, left=150, right=150):
    tcPr = cell._element.get_or_add_tcPr()
    tcMar = OxmlElement('w:tcMar')
    for m, val in [('top', top), ('bottom', bottom), ('left', left), ('right', right)]:
        node = OxmlElement(f'w:{m}')
        node.set(qn('w:w'), str(val))
        node.set(qn('w:type'), 'dxa')
        tcMar.append(node)
    tcPr.append(tcMar)

def set_callout_borders(cell, color_hex="1D63B8"):
    tcPr = cell._element.get_or_add_tcPr()
    tcBorders = OxmlElement('w:tcBorders')

    # Left thick border
    left = OxmlElement('w:left')
    left.set(qn('w:val'), 'single')
    left.set(qn('w:sz'), '24') # 3pt
    left.set(qn('w:space'), '0')
    left.set(qn('w:color'), color_hex)
    tcBorders.append(left)

    # Top, right, bottom nil
    for b_name in ['top', 'bottom', 'right']:
        b = OxmlElement(f'w:{b_name}')
        b.set(qn('w:val'), 'nil')
        tcBorders.append(b)

    tcPr.append(tcBorders)

def set_table_borders(table, top_sz=12, bottom_sz=12, header_sz=8):
    tblPr = table._element.xpath('w:tblPr')
    if not tblPr:
        return
    borders = OxmlElement('w:tblBorders')

    top = OxmlElement('w:top')
    top.set(qn('w:val'), 'single')
    top.set(qn('w:sz'), str(top_sz))
    top.set(qn('w:space'), '0')
    top.set(qn('w:color'), '000000')
    borders.append(top)

    bottom = OxmlElement('w:bottom')
    bottom.set(qn('w:val'), 'single')
    bottom.set(qn('w:sz'), str(bottom_sz))
    bottom.set(qn('w:space'), '0')
    bottom.set(qn('w:color'), '000000')
    borders.append(bottom)

    for b_name in ['left', 'right', 'insideH', 'insideV']:
        b = OxmlElement(f'w:{b_name}')
        b.set(qn('w:val'), 'nil')
        borders.append(b)

    tblPr[0].append(borders)

def build_docx(output_path):
    doc = docx.Document()

    # Set standard margins (1 inch = 1440 dxa)
    sections = doc.sections
    for section in sections:
        section.top_margin = Inches(1.0)
        section.bottom_margin = Inches(1.0)
        section.left_margin = Inches(1.0)
        section.right_margin = Inches(1.0)

    # Styles setup
    style_normal = doc.styles['Normal']
    style_normal.font.name = 'Arial'
    style_normal.font.size = Pt(10.5)
    style_normal.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
    style_normal.paragraph_format.line_spacing = 1.25
    style_normal.paragraph_format.space_after = Pt(6)

    # Helper functions for adding styled text
    def add_title(text):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(12)
        p.paragraph_format.space_after = Pt(12)
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(18)
        run.font.bold = True
        run.font.color.rgb = RGBColor(0x11, 0x3F, 0x7A)
        return p

    def add_subtitle(text):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_after = Pt(18)
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(11)
        run.font.italic = True
        run.font.color.rgb = RGBColor(0x55, 0x55, 0x55)
        return p

    def add_meta_info(text):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.LEFT
        p.paragraph_format.space_after = Pt(12)
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(9.5)
        run.font.color.rgb = RGBColor(0x44, 0x44, 0x44)
        return p

    def add_h1(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(18)
        p.paragraph_format.space_after = Pt(8)
        p.paragraph_format.keep_with_next = True
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(14)
        run.font.bold = True
        run.font.color.rgb = RGBColor(0x11, 0x3F, 0x7A)
        return p

    def add_h2(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(14)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.keep_with_next = True
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(12)
        run.font.bold = True
        run.font.color.rgb = RGBColor(0x00, 0x66, 0x66)
        return p

    def add_h3(text):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(10)
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.keep_with_next = True
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(11)
        run.font.bold = True
        run.font.color.rgb = RGBColor(0x33, 0x33, 0x33)
        return p

    def add_p(text, bold_prefix="", italic=False):
        p = doc.add_paragraph()
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.25
        if bold_prefix:
            run_b = p.add_run(bold_prefix)
            run_b.font.name = 'Arial'
            run_b.font.size = Pt(10.5)
            run_b.font.bold = True
            run_b.font.color.rgb = RGBColor(0x11, 0x3F, 0x7A)
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(10.5)
        run.font.italic = italic
        run.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
        return p

    def add_bullet(text, bold_prefix=""):
        p = doc.add_paragraph(style='List Bullet')
        p.paragraph_format.space_after = Pt(4)
        p.paragraph_format.line_spacing = 1.2
        if bold_prefix:
            run_b = p.add_run(bold_prefix)
            run_b.font.name = 'Arial'
            run_b.font.size = Pt(10)
            run_b.font.bold = True
            run_b.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
        run = p.add_run(text)
        run.font.name = 'Arial'
        run.font.size = Pt(10)
        run.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
        return p

    def add_callout(title, paragraphs):
        table = doc.add_table(rows=1, cols=1)
        table.alignment = WD_TABLE_ALIGNMENT.CENTER
        cell = table.cell(0, 0)
        cell.width = Inches(6.5)

        set_cell_background(cell, "F0F4F8") # Very light blue-grey
        set_cell_margins(cell, top=140, bottom=140, left=200, right=150)
        set_callout_borders(cell, "113F7A")

        p = cell.paragraphs[0]
        p.paragraph_format.space_after = Pt(6)
        r_t = p.add_run(title)
        r_t.font.name = 'Arial'
        r_t.font.size = Pt(11)
        r_t.font.bold = True
        r_t.font.color.rgb = RGBColor(0x11, 0x3F, 0x7A)

        for item in paragraphs:
            p_item = cell.add_paragraph()
            p_item.paragraph_format.space_after = Pt(4)
            p_item.paragraph_format.line_spacing = 1.2
            if isinstance(item, tuple):
                b_prefix, txt = item
                r_b = p_item.add_run(b_prefix)
                r_b.font.name = 'Arial'
                r_b.font.size = Pt(10)
                r_b.font.bold = True
                r_b.font.color.rgb = RGBColor(0x22, 0x22, 0x22)
                r_txt = p_item.add_run(txt)
                r_txt.font.name = 'Arial'
                r_txt.font.size = Pt(10)
            else:
                r_txt = p_item.add_run(str(item))
                r_txt.font.name = 'Arial'
                r_txt.font.size = Pt(10)

        # Add empty spacing paragraph after callout table
        p_space = doc.add_paragraph()
        p_space.paragraph_format.space_after = Pt(6)

    # --- DOCUMENT CONTENT GENERATION ---

    add_title("“衰老相关”在于标签，而非实测：\n肾脏衰老读数基因集的无泄露基准评估与临床验证")
    add_subtitle("‘Ageing-related’ by label, not by test: a leakage-free benchmark and clinical evaluation of gene sets used as kidney ageing readouts")

    add_meta_info("作者：刘孟晗¹,*，罗雅馨²,³，刘白客⁴，李贵森¹,*\n"
                  "单位：¹ 电子科技大学医学院、四川省医学科学院·四川省人民医院肾脏内科；² 英国布里斯托大学；³ MRC 综合流行病学研究所；⁴ 瑞士苏黎世大学\n"
                  "通讯作者：李贵森 (guisenli@163.com); 刘孟晗 (lmhdoc@163.com)")

    # 1. 摘要
    add_h1("一、 中英文摘要（Abstract）")

    add_h2("中文摘要")
    add_p("被标记为“衰老相关”（ageing-related）或“细胞衰老相关”（senescence-related）的基因集在肾脏病学研究中被广泛使用，用于将慢性肾脏病（CKD）描述为加速衰老、评估供体器官以及筛选抗衰老/清除衰老细胞（senolytic）靶点。然而在大多数研究中，这一标签仅仅记录了基因的来源（即来自某个衰老或细胞衰老数据库），而非它们在人类肾脏组织中随年龄变化的证据。当我们最初利用标准荟萃分析（meta-analysis）筛选策略从 9 个公开队列构建肾脏衰老签名时，它在发现集里表现强劲，但一旦核对样本记录并在保留（held-out）肾脏组织中进行验证，结果便不再成立。这一反差促使我们思考：整个领域内已发表的相关签名表现是否同样存在类似问题？", bold_prefix="背景：")

    add_p("从我们自建的 94 基因荟萃分析签名出发，我们首先评估了样本入组资格、循环评估（circular evaluation）、尸体解剖组织以及基因间相关性对表面衰老信号的影响。随后，我们在未参与签名构建的正常肾脏转录组中，对 48 个已发表的基因集（12 个源自人类肾脏年龄、7 个通用衰老资源、29 个在肾脏病研究中被标记为衰老/细胞衰老相关）进行了严格评估：包括 6 个供体队列（441 例供体）、TCGA 正常组织（n = 129）以及 CPTAC 瘤旁相邻组织（n = 75），GTEx 肾皮质（n = 85）则单独分析。在预先哈希存证（hashed plan）的分析计划下，我们将这些基因集与表达量匹配的随机基因进行比较（包含与不含瘤旁组织），并在基因级和评分级进行了组织成分矫正（基于 KPMP 单细胞核参考去卷积验证），同时开展了统计功效模拟以及 47 个 GTEx 组织的跨组织扩展。最后，我们在涵盖 2 个 CKD 队列（欧洲肾脏 cDNA 库 ERCB，n = 190，涵盖 9 种肾病；GSE175759，n = 87，含 eGFR）及 1 个移植前供肾活检队列（GSE147451，n = 151，含移植后肾功能结局）的共 428 例患者活检样本中直接验证了其临床效应。", bold_prefix="方法：")

    add_p("在我们最初的签名中，9 个公开数据集里有 3 个并非独立的正常肾皮质；若在保留队列之外重新进行基因筛选，年龄相关性几乎折半（ρ = 0.37 对比 0.20）。当我们转向已发表文献时，同样的鸿沟在更大范围内展现：作为 29 个疾病模型基因集中至少半数基因来源的 GenAge 和 CellAge 数据库，其在正常肾脏中的年龄相关性并不高于匹配的随机基因（富集比 0.94–1.04），且 CellAge 的诱导/抑制方向虽然在 47 个 GTEx 组织中的 29 个里追踪年龄，但在肾脏中却失效。29 个疾病模型基因集中有 5 个达到了 FDR < 0.05（剔除瘤旁组织后为 6 个），但因其中位基因数仅为 5 个，该检验对 30% 的富集效应仅有 25% 的统计功效。当把所有 29 个基因集汇总（共 151 个基因）时，其微弱的年龄富集（1.16–1.20）在基因级成分矫正后完全消失（富集比 1.05，P = 0.19；FDR < 0.05 的基因集为 0 个）；相比之下，源自肾脏年龄的基因集则保留了独立于成分的富集（汇总富集比 1.09–1.11，P = 0.0005；与 CKD 一致的肾脏年龄基因富集比达 1.61）。有符号（signed）肾脏年龄签名成功复现（ρ = 0.16–0.30），其中 45–57% 的信号归因于组织成分变化（与单细胞核 MuSiC 和 Bisque 去卷积高度一致，ρ = 0.87–0.95）。在患者活检样本中，GenAge、CellAge 及汇总疾病模型基因均在 CKD 对比参考组织中显著富集（富集比分别为 1.15、1.08 和 1.28；P = 0.0005），衰老签名在 9 种 CKD 病因中紧密追踪近端肾小管损伤（ρ = 0.72–0.89），且在校正肾小管损伤和组织成分后，消除了 CKD 活检样本中 70–100%（中位数 80%）的表面“加速衰老”效应。在移植前供体活检中，追踪年龄的签名对移植后肾功能均未提供超出供体实际日历年龄之外的额外预测价值（|ρ_adj| ≤ 0.05），而炎症与细胞衰老评分在校正年龄和成分后均失去了显著性。", bold_prefix="结果：")

    add_p("在肾脏病研究中，“衰老相关”仅描述了基因的数据库来源，而非其在肾脏中追踪年龄的证据。组织块（bulk）转录组肾脏年龄主要反映的是由成分改变驱动的肾小管损伤与炎症轴：它使 CKD 活检显得“更老”是因为组织受损，而在供体活检中则几乎无法提供超出实际年龄之外的新信息。R 包 kidneyAgeBench（v0.2.0）可在一个基因集投入使用前自动完成全部 5 项基准核验。", bold_prefix="结论：")

    add_p("肾脏衰老；基因签名；基准评估；统计功效；数据泄露；细胞衰老；单细胞核去卷积；慢性肾脏病；肾移植；转录组学", bold_prefix="关键词：")

    add_h2("English Abstract")
    add_p("Gene sets labelled \"ageing-related\" or \"senescence-related\" are widely used in kidney disease research to describe chronic kidney disease (CKD) as accelerated ageing, assess donor organs and nominate senolytic targets. In most studies, the label records where the genes came from (an ageing or senescence database) rather than evidence that they change with age in human kidney tissue. When we first built a kidney ageing signature from nine public cohorts using standard meta-analytic filters, it looked strong in the discovery series but did not hold up once we checked the sample records and tested it in held-out kidney tissue. That discrepancy prompted us to ask whether published signatures across the field behave any differently.", bold_prefix="Background: ")
    add_p("Starting from our 94-gene meta-analytic signature, we examined how sample eligibility, circular evaluation, post-mortem tissue and gene–gene correlation affected the apparent age signal. We then scored 48 published gene sets (12 derived from age in human kidney, 7 general ageing resources, 29 labelled ageing- or senescence-related in kidney disease studies) only in normal kidney transcriptomes not used to derive them: six donor cohorts (441 donors), TCGA normal (n = 129) and CPTAC tumour-adjacent tissue (n = 75), with GTEx cortex (n = 85) analysed separately. Sets were compared with expression-matched random genes under a hashed plan, before and after adjustment for tissue composition (validated against KPMP single-nucleus deconvolution), with power simulations and a 47-tissue GTEx extension. Clinical consequences were tested in 428 patient biopsies across two CKD cohorts (ERCB, n = 190; GSE175759, n = 87) and a pre-implantation transplant biopsy cohort (GSE147451, n = 151).", bold_prefix="Methods: ")
    add_p("In our initial signature, three of nine public series were not independent normal kidney cortex, and repeating gene selection outside each held-out cohort cut the age correlation nearly in half (ρ = 0.37 versus 0.20). GenAge and CellAge, the source of at least half the genes in 21 of 29 disease-model sets, were no more age-associated than matched genes in normal kidney (enrichment ratios 0.94–1.04), and CellAge's direction tracked age in 29 of 47 GTEx tissues but not in kidney. Five of 29 disease-model sets reached FDR < 0.05, but with a median of five genes the test had 25% power for a 30% enrichment. Pooled across all 29 sets (151 genes), their modest age enrichment (1.16–1.20) disappeared after gene-level composition adjustment (1.05, P = 0.19), whereas kidney age-derived sets retained composition-independent enrichment (pooled 1.09–1.11, P = 0.0005; CKD-concordant renal age genes 1.61). Signed kidney age signatures replicated (ρ = 0.16–0.30), with 45–57% attributable to composition. In patient biopsies, GenAge, CellAge and pooled disease-model genes were enriched for CKD versus reference tissue (ratios 1.15, 1.08 and 1.28; P = 0.0005), ageing signatures tracked proximal tubular injury across nine CKD etiologies (ρ = 0.72–0.89), and adjusting for tubular injury and composition removed 70–100% (median 80%) of the apparent \"accelerated ageing\" of CKD biopsies. In donor biopsies, age-tracking signatures added no prognostic value beyond donor chronological age (|ρ_adj| ≤ 0.05).", bold_prefix="Results: ")
    add_p("In kidney disease research, \"ageing-related\" describes where genes came from, not evidence that they track kidney age. Bulk transcriptomic kidney age largely reports a composition-driven tubular injury and inflammation axis: it makes CKD biopsies look older because they are injured and adds little beyond chronological age in donor biopsies. The kidneyAgeBench R package (v0.2.0) performs all five benchmark checks before a gene set is used.", bold_prefix="Conclusions: ")

    # 2. 背景
    add_h1("二、 背景（Background）")
    add_p("肾脏衰老发生早且肉眼可见。从第四个十年起，肾单位数量开始下降，肾硬化与间质纤维化不断累积，急性损伤后的修复能力持续减退 [1,2]。这些改变至少在三个临床场景中举足轻重：")
    add_bullet("在移植中，供体年龄是决定移植物存活率的最强因素之一，因此有研究提出通过供肾活检的转录组评分来提供超出实际日历年龄之外的信息 [3–5,15]。", bold_prefix="1. 肾移植评估：")
    add_bullet("在慢性肾脏病中，活检转录组被广泛用于论证 CKD 代表了一种早衰或加速衰老状态 [29]。", bold_prefix="2. 慢性肾脏病（CKD）：")
    add_bullet("在药物开发中，细胞衰老已成为治疗靶点，抗衰老/清除衰老细胞（senolytic/senomorphic）化合物的临床试验需要能够真实反映衰老细胞负荷是否改变的读数 [19]。", bold_prefix="3. 靶向抗衰老疗法：")
    add_p("上述每种应用都需要一个对肾脏衰老的精准测量：即在人类肾脏组织中随年龄变化，且能够与急性或慢性组织损伤相区分的指标。")
    add_p("公开的转录组数据大量提供了此类候选测量指标，其中一种构建范式占据统治地位：在糖尿病肾病（DKD）、狼疮性肾炎、IgA 肾病（IgAN）、移植排斥反应或肾细胞癌等队列中，将疾病组织与参考组织之间的差异基因与 GenAge 或 CellAge 等通用衰老/细胞衰老数据库求交集 [7]，随后通过机器学习对交集基因进行排序，提炼出“衰老相关”或“细胞衰老相关”的枢纽基因（hub genes）。所得基因列表随后被用作衰老读数：用于论证细胞衰老驱动疾病进展、构建预后模型或提名候选药物。")
    add_p("这一标签包含了一个隐式主张——即选出的基因会在人类肾脏组织中随年龄改变，但该标签实际上却是根据“出身”（provenance）赋予的。仅仅因为属于某个衰老数据库，就被直接等同于衰老生物学的证据。在 GTEx 跨组织数据中，根据 CellAge 诱导细胞衰老的基因确实倾向于随年龄表达更高 [25]，基于表达量构建的细胞衰老签名也与衰老同向改变 [26]；细胞衰老基因集也被用作评估衰老细胞的表型标记 [27]。然而，肾脏病研究从这些数据库中衍生出的规模小得多的基因集是否能在肾脏组织中追踪年龄，此前从未得到检验。检验这一点比表面上看起来困难得多：")
    add_bullet("在产生数据的同一批数据中打分，会让基因集看起来比实际更强 [8]；", bold_prefix="• 循环评估偏性：")
    add_bullet("显著的相关性并不代表特异性，因为任意随机基因集也常与临床变量相关，因此受检基因集必须与表达量匹配的随机基因相比（而非与全基因组相比）[9]，并需针对基因间相关性进行评估 [30]；", bold_prefix="• 缺乏特异性对照：")
    add_bullet("衰老会改变肾皮质的细胞组成，因此一个基因集可能仅仅因为免疫、基质或肾小管含量的比例漂移而表现出“追踪年龄”[31]；", bold_prefix="• 细胞成分混杂：")
    add_bullet("“正常”肾脏的公开来源（活体供体活检、切缘正常组织及尸体解剖皮质）各自带有特定的偏性；", bold_prefix="• 组织来源异质性：")
    add_bullet("大多数疾病衍生基因集仅包含极少数基因，小基因集未达到统计显著并不必然代表其与年龄无关。", bold_prefix="• 统计功效不足：")
    add_p("我们正是因为尝试构建肾脏衰老签名而关注到这一问题的。利用 9 个公开正常肾脏队列的标准荟萃分析流程，我们识别出 94 个基因，其复合评分在发现集里与供体年龄展现出良好的相关性（ρ = 0.37）。然而，当我们仔细核对底层的样本记录并在挑选基因之外的队列中进行验证时，结果却不复成立：3 个公共数据集并非独立的正常肾皮质，交叉队列保留评估使效应几乎折半，评分在 GTEx 尸体皮质中完全平坦，而看似强劲的基因级一致性信号在控制基因间共表达后瞬间消失。如果直接从供体年龄构建的签名都会出现这种假象，那么肾脏病研究中仅仅因为与数据库重叠就被冠以“衰老相关”的众多已发表基因集又当如何？为了找出答案，我们在未参与其构建的正常肾脏转录组中对 48 个已发表基因集进行了评分，在预先哈希存证的分析计划下将每个基因集与表达量匹配的随机基因进行比较，估计了在该基因集规模下的统计功效，在评分级和单基因级将年龄相关性与组织成分分离开来（基于单细胞核参考去卷积验证 [31–33]），评估了衰老数据库在 47 个其他人类组织中是否表现不同，并在涵盖 9 种 CKD 病因和肾移植的 428 例肾活检样本中直接测试了其临床后果。我们将该基准评估及全部 5 项使用前核验打包为 R 包 kidneyAgeBench 发布。")

    # 3. 方法
    add_h1("三、 方法（Methods）")

    add_h2("3.1 研究设计与分析计划")
    add_p("本研究包含四个部分：(1) 在公开队列中构建并测试荟萃分析肾脏衰老签名；(2) 在正常人类肾脏中对 48 个已发表基因集进行冻结基准评估（locked benchmark）；(3) 事后（post hoc）方法学扩展，解析统计功效、组织成分及组织特异性；(4) 在 CKD 及移植活检样本中进行直接临床验证（图 1a）。基准评估计划、基因集库及打分代码在任何库内基因集被打分前均已冻结，并附加 SHA-256 哈希时间戳。计划明确了评估数据集、数据泄露规则、评估指标以及三项预期：(E1) 疾病模型基因集的中位富集比低于 1.2，且达到 FDR < 0.05 的比例不足一半；(E2) 源自肾脏年龄的有符号基因集的复现率高于疾病模型基因集的“衰老信息度”；(E3) 组织成分校正会将有符号相关性削弱中位数 30% 以上。分析计划预先披露了部分基因集在部分数据集中的早期打分历史。在初始运行后，使用完全相同的代码追加了两次基因库扩展（分别为 9 个和 5 个基因集），且每次扩展均在打分前完成哈希存证。初始冻结运行之后追加的分析均标注为事后分析（post hoc）。")

    add_h2("3.2 评估数据集")
    add_p("我们最初收集了 9 个包含供体年龄的 GEO 正常或参考人类肾脏数据集。每个数据集均从原始数据重新构建，样本元数据对照原始文献重新核对。入组规则为：具备全基因组表达谱的独立成年正常/参考肾皮质供体。3 个数据集未满足此规则（表 1）：GSE232249 包含来自 6 颗肾脏的 m6A 免疫沉淀和上清 RNA 片段；GSE21785 包含来自 6 例活体供体的肾小管和肾小球微切割样本；GSE217427 包含来自 16 例供体的 22 例皮质样本（8 例伴急性肾损伤，且存在双侧采样）。在 GSE147451 中，177 例移植前活检来自 166 例遗体供体，11 例贡献双肾的供体被平均合并为每例供体一个表达谱。满足条件的 6 个队列（GSE36895、GSE60807、GSE147451、GSE66272、GSE93069 和 GSE167093）共计包含 441 例供体，统称为 G6。", bold_prefix="GEO 供体队列：")
    add_p("TCGA 肾透明细胞癌（KIRC）、肾乳头状细胞癌（KIRP）和肾嫌色细胞癌（KICH）项目的癌旁正常组织（Solid Tissue Normal）从 GDC 下载为 STAR 基因级 TPM。每例患者保留 1 个样本（n = 129；KIRC 72 例、KIRP 32 例、KICH 25 例；年龄 28–89 岁；女性 42 例）。表达量分析采用 log2(TPM + 1)，保留在至少 20% 样本中表达的基因。", bold_prefix="TCGA 队列：")
    add_p("来自 CPTAC 肾透明细胞癌发现队列的癌旁相邻组织（NAT）RNA 测序表达谱（上四分位数归一化，log2(x + 1)）与泛癌临床表中的年龄和性别相链接（n = 75）。", bold_prefix="CPTAC 队列：")
    add_p("GTEx v8 肾皮质 TPM 表达谱（n = 85）以十年年龄段的中点进行编码分析。GTEx 单独汇报，因其为带有粗粒度年龄信息的尸体解剖组织。", bold_prefix="GTEx 队列：")
    add_p("为测试衰老与细胞衰老签名在临床活检中的表现，我们分析了 3 个具备临床终点的队列：(1) 欧洲肾脏 cDNA 库（ERCB）肾小管间质转录组队列（GSE104954；n = 190：21 例活体供体及 169 例涵盖 9 种 CKD 诊断的患者）；(2) 独立微切割肾小管间质 RNA 测序队列（GSE175759；n = 87：22 例对照及 65 例 CKD 活检，均带有 eGFR 记录）；(3) GSE147451 中带有移植后肾功能结局记录的 151 例遗体供体。", bold_prefix="临床 CKD 与移植队列：")

    add_h2("3.3 构建与测试标准荟萃分析肾脏衰老签名")
    add_p("在最初收集的 9 个数据集各自内部，估计每个基因与年龄的偏相关（校正性别，如可用）并转化为 Fisher z 值。基因级估计值通过限制最大似然（REML）随机效应荟萃分析（metafor）进行合并。一个基因若在至少 8 个数据集中测得、达到 FDR < 0.05、在至少 70% 的数据集中方向与汇总方向一致，且在依次剔除每个数据集时仍保持 P < 0.05 的同向性，则进入签名。这产生了 94 个基因（KAS-94）。将同样规则应用于 G6 产生 82 个基因（KAS-82）。为了量化选择乐观偏性，我们将签名评分与年龄的样本内合并相关性，与在保留队列之外重复基因筛选的嵌套留一队列估计值进行了对比。在 GTEx 中，我们计算了 KAS-94 评分与年龄的相关性，并以二项式检验和置换检验测试了方向一致性。")

    add_h2("3.4 基因集库与数据泄露控制")
    add_p("基因集被分为三类（附表 S1）：A 类（n = 12）源自人类肾脏组织年龄或供体质量；B 类（n = 7）为通用衰老或细胞衰老资源（如 GenAge, CellAge, SenMayo）；C 类（n = 29）为肾脏病研究筛选或标记为衰老/细胞衰老相关的基因集。数据泄露控制规则：任何数据集均不得用于评估由其衍生出的签名。因此，KAS-94 和 KAS-82 仅在 G6 之外进行评估，Rowland 签名仅在 TCGA 之外评估。")

    add_h2("3.5 基因级富集、成分校正与统计功效模拟")
    add_p("在每个数据集内部，估计每个基因与年龄的偏 Spearman 相关。Fisher z 统计量通过 Stouffer 方法在无泄露的数据集中合并。基因集的统计量为其测得成员的平均绝对合并 z 值。零分布来自 2,000 个表达量均值和标准差十分位数匹配的随机基因集。我们汇报富集比（观察统计量除以零分布均值）、经验 P 值以及 BH FDR。协主分析从年龄统计量和零分布中均剔除了瘤旁相邻组织。为分离成分影响，我们在校正 5 个标志物组织成分评分（免疫、基质、内皮、远端小管和近端小管含量）后重新计算了基因级年龄相关性。同时，我们在基因集规模（2–100 个基因）下模拟了基准检验的统计功效，导出了 80% 功效下可检测到的最小富集比（mde80）。")

    add_h2("3.6 临床评估与 GTEx 跨组织扩展")
    add_p("在 ERCB (n=190) 和 GSE175759 (n=87) 中，估计了 CKD 对比参考活检的表面衰老偏移，并在校正组织成分和 10 基因近端小管损伤评分前后评估了相关性。在 GSE147451 (n=151) 中，测试了各签名与移植后肾功能不良的关联。此外，从 GTEx v8 获取 47 个组织的 TPM 文件，对 CellAge 诱导/抑制方向在跨组织中的追踪能力进行了检验。")

    # 4. 结果
    add_h1("四、 结果（Results）")

    add_h2("4.1 标准构建流程揭示的四个陷阱")
    add_p("在发现系列中，自建的 KAS-94 签名复合评分与供体年龄的相关性达到合并 ρ = 0.37（95% CI 0.29–0.45；图 2a）。但重新核对样本描述显示 3 个系列不符合独立正常肾皮质要求（表 1，图 1c）。清理后缩减为 G6 队列（441 例供体），重新筛选得到 KAS-82。在嵌套留一队列循环中从头重复基因选择并在保留队列中打分，合并相关性降至 ρ = 0.20（图 2a）。带到 GTEx 尸体肾皮质（n = 85）时，评分与年龄完全无关（ρ = 0.04；图 2b）。单基因方向一致性在假设独立时 P = 4.7×10⁻⁶，而置换检验为 P = 0.20（图 2c）。")

    add_h2("4.2 已发表“衰老相关”基因集的基准评估")
    add_p("图 3 展示了 48 个基因集的评估结果。冻结主分析中，12 个 A 类基因集中有 9 个达到 FDR < 0.05（中位富集比 1.23）；29 个 C 类基因集中仅 5 个达到 FDR < 0.05（中位富集比 1.05，14 个富集比 < 1）。剔除瘤旁组织的协主分析中，A 类 7/10、C 类 6/29 达到 FDR < 0.05（中位富集比分别为 1.31 和 1.14；图 7c）。使用 camera 方差膨胀因子（VIF）校正基因间相关性后，达到 FDR < 0.05 的 C 类集在冻结和协主分析中均减少至 4 个。")

    add_h2("4.3 小基因集的统计功效限制")
    add_p("C 类基因集的中位测得基因数仅为 5 个。统计功效模拟显示（图 4a），在 5 个基因时，对 1.3 的真实富集比仅有 25% 的功效，80% 功效下所需的最小富集比（mde80）达 1.83。没有任何 C 类基因集具备 80% 功效来检测 1.3 的富集比（功效范围 16%–51%），因此单独 C 类基因集的非显著结果不具备明确的否定信息量。")

    add_h2("4.4 数据库标签不追踪肾年龄，成分校正消除 C 类基因信号")
    add_p("在汇总基因源层面（功效 >90%），GenAge 和 CellAge 在正常肾脏中与年龄的相关性并不高于匹配随机基因（富集比 0.94 和 0.99；图 8a）。21 个 C 类基因集至少半数基因来自这两个数据库。成分校正前，151 个汇总 C 类基因表现出微弱的年龄相关（富集比 1.16）；但在校正 5 个区室标志物评分后，其富集完全消失（富集比 1.05, P = 0.19；0/29 个 C 类集 FDR < 0.05）。相反，外部 A 类基因池保留了显著的成分独立富集（富集比 1.09, P = 0.0005；CKD 一致的 31 个肾年龄基因富集比达 1.61）。跨组织扩展显示（附图 S2），CellAge 方向在 29/47 个 GTEx 组织中追踪年龄，但在肾脏中失效（超额 z = 0.87, P = 0.20）。")

    add_h2("4.5 表达量梯度与组织成分漂移造成假象")
    add_p("与全基因组相比，C 类基因看似强相关（中位 |z| = 2.03 对比全基因组 1.20；图 5a）。这是因为高表达和高变异基因具有更高的年龄相关性，而 72% 的 C 类基因落在前 3 个表达十分位数中（图 5b）。成分校正平坦化了这一梯度，使全基因组中位 |z| 降至 0.81。跨 C 类模型最常重复使用的基因（如 FOS, CDKN1A, STAT1, EGF, PCK1；图 5d）多为随损伤与炎症改变的受损/应激基因。")

    add_h2("4.6 有符号签名的复现性与组织来源影响")
    add_p("除 MHC II 之外的有符号 A 类签名均显著超越规模匹配的随机签名（合并 ρ 0.16–0.30；图 6a）。成分校正将合并相关性削弱了中位数 57%（图 6c），且成分评分与 KPMP 单细胞核去卷积（MuSiC, Bisque）高度一致（ρ = 0.87–0.95；图 8b）。在组织来源方面，中位有符号 ρ 在 G6 为 0.18，TCGA 为 0.19，但在 CPTAC 瘤旁组织中高达 0.40（图 7b）。男性 Y 染色体基因评分在 CPTAC 中随年龄陡峭下降（ρ = -0.75），表明存在年龄依赖的肿瘤场效应（图 7a）。")

    add_h2("4.7 临床活检中的损伤报告与移植预测价值")
    add_p("在 ERCB (n=190) 和 GSE175759 (n=87) 中，GenAge、CellAge 和汇总 C 类基因在 CKD 对比参考肾脏中强富集（富集比 1.15, 1.08 和 1.28, P = 0.0005；图 8a）。在 9 种 CKD 病因中，所有衰老签名均与近端小管损伤评分紧密耦合（ρ = 0.72–0.89；图 8c）。校正成分和小管损伤消除了 CKD 活检样本中 70–100%（中位数 80%）的表观“加速衰老”偏移（图 8d）。在 151 例移植前遗体供体活检中，供体实际日历年龄与移植后不良 GFR 相关的 ρ = 0.30；校正供体日历年龄和成分后，所有衰老和细胞衰老签名对移植后结局均无额外预测价值（|ρ_adj| ≤ 0.05；图 8e）。")

    # 5. 讨论与结论
    add_h1("五、 讨论与结论（Discussion & Conclusions）")

    add_h2("5.1 讨论")
    add_p("对于肾脏病研究中被冠以“衰老相关”的大多数基因集，该标签仅是对“出身”的陈述：即这些基因来自某个衰老或细胞衰老数据库，而非它们在肾脏中随年龄改变的证据。单独的标签基因集规模过小（中位数 5 个基因），使其非显著结果本身缺乏否定性结论的明确信息量；而当把 29 个基因集的 151 个基因合并时（功效 >90%），其微弱的未校正年龄相关在校正组织成分后彻底消失（富集比 1.05, P = 0.19）。相比之下，衍生自肾脏年龄的签名保留了独立于成分的年龄信号，特别是与 CKD 一致的基因（富集比 1.61）。")
    add_p("三条观察解释了疾病基因是如何获得“衰老基因”外貌的：(1) 表达量水平与组织成分耦合，导致高表达区室转录本在组织块中展现出表面年龄相关；(2) 数据库出身在疾病对比正常衰老中的混淆，即受损肾脏激活了通用应激与细胞周期基因；(3) 基因间相关性低估了朴素检验的零分布方差。有符号肾脏年龄签名的合并相关性为 ρ = 0.16–0.30，其中约一半反映了组织成分漂移。在组织块转录组中，人类肾脏衰老主要表现为一个低度慢性肾小管损伤和炎症轴，而非一套独特的长寿程序。")

    add_h2("5.2 临床应用启示")
    add_bullet("在 CKD 中升高的衰老评分反映的是小管间质损伤和白细胞浸润，而非独立加速的衰老钟；疾病队列中早衰肾脏的主张需要提供校正损伤和成分后的分析。", bold_prefix="1. 慢性肾脏病（CKD）：")
    add_bullet("在使用分子“年龄”评分来接受、分配或弃用供体肾脏之前，应当检验其是否具备超出供体日历年龄和活检成分之外的增量价值。", bold_prefix="2. 肾移植评估：")
    add_bullet("在临床试验中使用组织块衰老签名作为药效学终点，面临着测得的是抗炎或成分改变而非衰老细胞清除的风险，需要单细胞或空间分辨读数。", bold_prefix="3. 抗衰老疗法临床试验：")

    add_h2("5.3 局限性与结论")
    add_p("本研究的局限性包括：正常评估数据为组织块转录组；部分扩展分析为事后设计；功效模拟基于特定加权假定；大部分 C 类模型最初训练用于疾病分类而非衰老预测；以及缺乏细胞内在生物学肾脏年龄的独立金标准。")
    add_p("结论：基因集“衰老相关”标签仅代表其数据库来源。组织块转录组衰老读数主要测量肾小管损伤与炎症轴，使 CKD 活检显得更老是因为组织受损，而在供体活检中几乎无法提供超出日历年龄之外的信息。在使用基因集作为肾脏衰老读数前，应通过 kidneyAgeBench 完成全部 5 项核验。", bold_prefix="结论：")

    # 6. Box 1与图表说明
    add_h1("六、 关键提炼框与图表说明（Box 1 & Figure Legends）")

    # Callout for Box 1
    box1_items = [
        ("1. 无泄露肾脏年龄检验：", "证明基因在未参与挑选它们的正常人类肾脏组织中随年龄改变。"),
        ("2. 感知表达量与相关性的零分布：", "与表达量匹配的随机基因集（并控制基因间相关性）进行比较，而非全基因组或朴素二项式检验。"),
        ("3. 成分与损伤分离：", "在汇报未校正估计值的同时汇报成分校正估计值（并在 CKD 活检中校正肾小管损伤）。"),
        ("4. 组织来源健壮性：", "声明每个评估数据集的组织来源，并同时汇报包含与剔除瘤旁或尸体解剖组织时的估计值。"),
        ("5. 在汇报基因集规模下的统计功效：", "汇报测得基因数、检测现实富集（如 1.3）的统计功效以及 mde80；将小基因集的非显著结果视为“无明确信息量”。"),
        ("", "在 kidneyAgeBench (v0.2.0) 中，单次调用 kab_benchmark(genes, direction, exclude, config) 即可自动完成全部 5 项核验。")
    ]
    add_callout("Box 1. 将一个基因集描述为肾脏衰老读数之前的五项核验", box1_items)

    add_h2("图表说明（Figure Legends）")
    add_p("workflow。48 个已发表基因集在正常肾转录组中打分（剔除衍生数据集），随后在 428 例临床活检中评估；分析计划、基因库及代码在打分前均已哈希存证。(b) 正常肾评估数据集中的供体年龄。点代表供体，横线代表中位数，标签给出供体数及女性百分比。(c) 最初收集的 9 个 GEO 系列的入组资格审查（表 1）。", bold_prefix="图 1. 研究设计与评估数据。(a) ")
    add_p("KAS-94 和 KAS-82 与年龄的合并偏 Spearman ρ (95% CI)：选择队列中的样本内估计值（空心圆）与嵌套留一队列估计值（实心圆）。(b) 未参与构建的 GTEx 肾皮质 (n = 85) 中的 KAS-94 评分与年龄关系。(c) 93 个测得的 KAS-94 基因中 GTEx 年龄相关性具备发现方向的基因数，与假设独立的二项式零分布（红）及置换检验零分布（灰）相比。", bold_prefix="图 2. 常规构建的肾脏衰老签名在训练流程之外测试表现。(a) ")
    add_p("48 个基因集的评估计分卡（测得基因数、基因级富集比、1.3 富集比下的统计功效、标称 P < 0.05 比例、方向一致性比例、有符号合并 ρ 及成分校正前后相关性）。(b) 在冻结主分析和协主分析中，各类别达到 FDR < 0.05、高于随机但未达显著、或低于随机的基因集比例。", bold_prefix="图 3. 48 个肾脏衰老、衰老与细胞衰老基因集的基准评估计分卡。(a) ")
    add_p("冻结配置下按测得基因数的模拟统计功效（经验 P < 0.05，真实富集比 1.1–2.0）。(b) 各 C 类基因集的富集比对其规模；虚线为 80% 功效下最小可检测富集比。(c) 汇总基因源相对于表达量匹配基因的平均绝对年龄 z（95% 零分布区间）。(d) 如果所有 29 个基因集共享相同真实富集比时预期的 P < 0.05 数量对比观察值。", bold_prefix="图 4. 基因级检验的统计功效与汇总证据 (post hoc)。(a) ")
    add_p("所有测得基因及各基因池的绝对合并年龄 z；虚线标记全基因组中位数。(b) CPTAC 癌旁 RNA 各表达量均值和标准差十分位数分箱的平均绝对年龄 z；圆圈显示 C 类基因数。(c) 各 C 类基因集中同时属于 GenAge 或 CellAge 的基因百分比对其富集比。(d) 在 2 个或更多 C 类模型中使用的基因的合并年龄 z；灰色带标记 |z| < 1.96。", bold_prefix="图 5. 为何“衰老相关”基因集看起来与年龄相关 (post hoc)。(a) ")
    add_p("跨无泄露数据集的签名评分与年龄合并偏 Spearman ρ (95% CI，菱形) 及 500 个规模匹配随机签名的分布（灰色）。(b) 每个队列中的偏 Spearman ρ；点标记 P < 0.05。(c) 基于标志物的组织成分校正前后的合并 ρ。(d) 每个数据集中每个成分评分与年龄的偏 Spearman ρ。", bold_prefix="图 6. 有符号签名的复现与组织成分。(a) ")
    add_p("男性中 9 个 Y 染色体基因评分与 TCGA 正常组织、GTEx 皮质及 CPTAC 瘤旁组织中年龄的关系。(b) 按组织来源划分的有符号评分与年龄的 ρ。(c) 冻结主分析中每个基因集的富集比对其在协主分析中的富集比。", bold_prefix="图 7. 组织来源塑造表面衰老信号。(a) ")
    add_p("跨三项评估的汇总基因源富集（未校正正常肾年龄、成分校正正常肾年龄、以及 CKD 对比参考活检 n = 277）。(b) G6 年龄相关性（n = 414）在标志物校正后与 KPMP 单细胞核参考去卷积（MuSiC, Bisque）的一致性。(c) ERCB (n = 190) 中 9 种 CKD 病因的中位签名评分对中位近端小管损伤评分作图。(d) ERCB CKD 活检对比活体供体的表面衰老偏移在校正成分和损伤前后对比。(e) 151 例移植前供体活检中签名与移植后肾功能不良的 Spearman 相关（校正年龄及成分前后）。", bold_prefix="图 8. 基因级成分解析、单细胞核去卷积一致性及临床评估 (post hoc)。(a) ")

    # Save document
    doc.save(output_path)
    print(f"Successfully generated {output_path}")

if __name__ == "__main__":
    import sys
    out_file = sys.argv[1] if len(sys.argv) > 1 else "manuscript_v6_GenomeMed_ZH.docx"
    build_docx(out_file)
