# Daniel 2022 (GSE188666): the declaration

Phase 3 Task 6, 2026-11-20. Declared before the derived Daniel file is built and before any method
runs on it. Written from three things only: Task 3's note with Task 5a's sections
([2026-11-19-daniel2022-metadata.md](2026-11-19-daniel2022-metadata.md)), the GSM records that note
transcribes, and the literature below. No contig file and no metadata file was opened for this task;
the label counts used for the "more than 50 cells" rule are the per-capture `ident` counts that
Task 3 recorded (note section 3).

**No clone was counted and no method was run on Daniel 2022 before this declaration.** No clone x
state, clone x capture or clone x lineage count was computed, and no clone count of any kind.

The decision rules are the plan's
([phase 3 plan](../plans/2026-11-12-phase-3-tcell-data.md), Task 6), applied in its order. The spec
(§6.1) wins over the plan where they differ; no conflict was found on the rules themselves (see
"For the controller" at the end).

## 1. Literature

Every DOI was checked against Crossref (`https://api.crossref.org/works/<doi>`) on 2026-11-20: title,
journal and year match for all six, so no correction was needed. Full texts were read as the NIH
author manuscripts in PMC, fetched as JATS XML through NCBI E-utilities
(`efetch.fcgi?db=pmc&id=<PMCID>`), which serves no CAPTCHA. Europe PMC's `fullTextXML` was tried
first and on retries returned HTTP 500 or 503, or timed out, for all five. No CAPTCHA or sign-in was
met or bypassed.

| Paper | DOI (Crossref, verified) | Claim used here | Method behind it | Where read |
|---|---|---|---|---|
| Hudson W.H. et al., Immunity 51:1043, 2019 | 10.1016/j.immuni.2019.11.002 | Stem-like Tcf-1+ CXCR5+ CD101−Tim3− cells give rise to a transitory CD101−Tim3+ population; the transitory cells convert into CD101+Tim3+ cells; CD101+Tim3+ cells stay as they are, even 28 days after transfer. The transitory cells preferentially express Cx3cr1, Klrg1, Tbx21, Zeb2, S1pr1, Klf2 and Mki67; CX3CR1 protein is found only on them. | Adoptive transfer of sorted subsets into infection-matched hosts (Results, "CD101−Tim3+ cells are a transitory CD8+ T cell population arising from Tcf-1+ stem-like cells"; Fig. 4); RNA-seq of sorted subsets (Fig. 3) | PMC6920571, full text (E-utilities) |
| Zander R. et al., Immunity 51:1028, 2019 | 10.1016/j.immuni.2019.10.009 | Three subsets at day 30: Ly108+ (Tcf7, Slamf6, Id3, Cxcr5; TCF-1+), CX3CR1+ (Klre1, Klra9, Klrd1, Klrg1, Tbx21, Zeb2, Klf2; KLRG1 and T-bet high by flow) and CX3CR1−Ly108− (PD-1, Tim3, 2B4, CXCR6, Eomes high). Ly108+ cells give rise to both other subsets. Most CX3CR1+ donor cells keep their phenotype and do not become Ly108+ or CX3CR1−Ly108−. About 40-50% of the progeny of CX3CR1−Ly108− donors acquire CX3CR1. | Adoptive transfer of sorted subsets into naive hosts then infected with Clone 13 (Results, "Slamf6+ CD8+ T cells give rise to both Pdcd1+ exhausted and Cx3cr1+ effector CD8+ T cell subsets"; Fig. 4). Its Monocle 2 trajectory is computational and is not used. | PMC6929322, full text (E-utilities) |
| Beltra J.-C. et al., Immunity 52:825, 2020 | 10.1016/j.immuni.2020.04.014 | Four stages: Texprog1 (Ly108+CD69+) and Texprog2 (Ly108+CD69−), both TCF1+ and PD-1 intermediate; Texint (Ly108−CD69−; Cx3cr1, Klrg1, Klr family, Gzma/b, Tbx21, Zeb2, Ki67); Texterm (Ly108−CD69+; Pdcd1, Lag3, Tigit, Cd244, Entpd1, Cd101, Cd38). Ly108+ donors regenerate all four; Ly108−CD69− donors mainly keep their identity and give rise to Ly108−CD69+; Ly108−CD69+ donors keep their state; neither Ly108− subset gives rise to Ly108+ cells. | Adoptive transfer of sorted, CFSE-labelled subsets into infection-matched hosts (Results, "Developmental Relationships between the Four Ly108- and CD69-Defined Tex Subsets"; Fig. 2); RNA-seq (Fig. 4) | PMC8360766, full text (E-utilities) |
| Im S.J. et al., Nature 537:417, 2016 | 10.1038/nature19330 | The PD-1+CXCR5+Tim-3− subset is TCF1+, self-renews and gives rise to the CXCR5−Tim-3+ subset; CXCR5−Tim-3+ cells keep their phenotype. | Adoptive transfer of sorted subsets into infection-matched hosts (Fig. 3) | PMC5297183, full text (E-utilities) |
| Utzschneider D.T. et al., Immunity 45:415, 2016 | 10.1016/j.immuni.2016.07.021 | A Tcf1-expressing subpopulation sustains the response in chronic infection and expresses PD-1 and Lag-3. | Abstract only | PubMed abstract (PMID 27533016, E-utilities). The full text is not in PMC or Europe PMC; Elsevier's text-mining API refused the request without a key. Its Results are **not read**, so only the abstract's claim is used, and only as support for the root. |
| Daniel B. et al., Nat Immunol 23:1614, 2022 | 10.1038/s41590-022-01337-5 | The marker description of each scRNA-seq cluster, the sorts behind the D21 spleen captures, and the Methods (section 2). Its clonal analysis (Figs. 3, 5, 6 and every "clone" statement) is not used for any mapping or lineage. | DEG-based cluster annotation (Fig. 1e, Supplementary Table 1), flow cytometry and adoptive transfers | PMC11225711, full text (E-utilities). Task 3 recorded it as unreadable through PMC's HTML (CAPTCHA) and Europe PMC (HTTP 500); the E-utilities API serves it. |

## 2. What each label is, from Daniel's non-clonal description

Daniel annotates 11 scRNA-seq clusters by differentially expressed genes (Results, "A multi-omic atlas
of T cell exhaustion"). The metadata's 11 `ident` labels (note section 3) match those clusters one to
one by name. `TEx` is Texterm: it is the one exhausted label without a qualifier, and it is
corroborated by Task 3's label counts (9,216 of 9,938 liver TetNeg cells are `TEx`; Daniel reports the
liver as 93.9% Texterm).

| `ident` label | Daniel's cluster | Daniel's markers (non-clonal) |
|---|---|---|
| Naive | Tnaive | Ccr7, Sell, Lef1 |
| Effector | Teff | Klrg1, Ly6c2 (mainly acute) |
| Effector_Memory | Tem | Klrb1c, Klrd1, S1pr1 (acute) |
| Memory | Tmem | Il7r, Arl4c, Il18r1 (acute) |
| TEx_prog | Texprog | Tcf7, Slamf6, Id3; sorted as PD-1+SLAMF6+ |
| TEx_int | Texint | Lgals3, S100a4, Mki67; Ctla4, Pdcd1, Lag3; intermediate CXCR6; by flow CX3CR1+SLAMF6−KLRG1− |
| TEx | Texterm | Gzma, Cd101, Entpd1; Cxcr6, Cd101, Tigit higher than in TexKLR; by flow CX3CR1−SLAMF6− |
| TEx_KLR | TexKLR | S1pr5, Cx3cr1, Klrc1; KLR family (Klrd1, Klrk1, Klrc1, Klre1, Klrg1) and Zeb2 specific to it; by flow CX3CR1+KLRG1+SLAMF6− |
| TEx_earlyEffector | Texeeff | Xcl1, Top2a, Mif; predominant at D8 |
| TEx_lung | Texlung | Lag3, Ifng, Ccl4; a lung variant of terminal exhaustion |
| Ifn_Response | TexISG | Isg15, Ifit1, Isg20 |

**The CX3CR1+ subset of the three papers is two Daniel labels.** Daniel opens the section titled
"CX3CR1+ Tex comprise two functionally distinct subsets" by citing Hudson and Zander (its refs 8 and 9)
for CX3CR1+ cells as an intermediate between Texprog and Texterm, then reports that sorted CX3CR1+
cells span Texint and TexKLR. Task 3's counts agree: the D21 spleen PD-1+CX3CR1+ sort (GSM5687807's
cells) holds 2,805 `TEx_KLR` and 1,521 `TEx_int` cells. The papers' markers do not pick one of the two
either: Hudson's transitory cells and Beltra's Texint carry Mki67/Ki67 (Daniel's Texint marker) and
also Klrg1, the Klr genes and Zeb2 (Daniel's TexKLR markers); Zander's CX3CR1+ subset carries the
KLR and Zeb2 signature, but the population whose fate Zander tested is the CX3CR1+ sort gate, which
by Daniel's flow definitions holds both TexKLR (KLRG1+) and Texint (KLRG1−) cells. So each paper's
CX3CR1+ subset maps to `TEx_int` and `TEx_KLR` together, not to one label.

## 3. The rules, in order

**Arm.** Chronic (LCMV-c13) captures only (spec §6.1). Candidates: the 15 chronic captures (liver and
lung D21: GSM5687816-19; spleen D21: GSM5687824-32; spleen D8: GSM5687833-34; spleen D5:
GSM6381666-69). The 4 acute captures are out.

**Days: 8 and 21; Day 5 excluded.** The main matrix holds chronic days 8 and 21. The D5 matrix is
kept only if the records say both matrices were processed together *and* its labels share the main
vocabulary. Neither holds. The records say nothing on joint processing (note section 6: the D5 GSMs
were submitted eight months later, name no processed file, and the D5 metadata is a separate Seurat
object with its own embedding), and Daniel's Methods list "Chronic including Day 5" as its own
clustering with its own parameters. The D5 labels are a different vocabulary of 7 (it adds
`TEx_earlyEffector2` and lacks five main labels), and on the 43,820 cells in both objects the labels
were assigned apart (note section 3). Reason: "a matrix processed apart would put a batch on the day
axis". So GSM6381666-69 are out, and the D5 matrix's copies of the D8/D21 cells are never read.

**Tissues: spleen only.** Liver and lung were sampled only at D21 (note section 3); spleen at D8 and
D21. Keeping liver or lung would confound tissue with day. GSM5687816-19 are out. Candidates left: the
11 chronic spleen D8/D21 captures, GSM5687824-34.

**State gates.** A sort on Slamf6/Ly108, CX3CR1, TIM-3, CD101, KLRG1, CXCR5 or TCF-1 is a state gate;
tetramer and PD-1 gates are not (note section 4).
- State-gated: GSM5687824 (CD8+ CX3CR1− SLAMF6+), GSM5687825 (CD8+ CX3CR1−), GSM5687826 (CD8+ CX3CR1+).
- Not state-gated: GSM5687827 (CD8+ PD1+; PD-1 gate), GSM5687828 (gp33 mid), GSM5687829 and GSM5687830
  (gp33−), GSM5687831 and GSM5687832 (gp33+), GSM5687833 (D8 gp33−), GSM5687834 (D8 gp33+).
  Daniel's Results add that the D21 spleen subset sorts were PD-1+ with SLAMF6/CX3CR1 gates; that does
  not change which sorts carry a state marker.

**Namespace: rule (b), namespace = capture, assessable by this rule.** (a) needs the records to name
each capture's mouse; they name no mouse (note section 5: no GSM or BioSample field), so (a) cannot
apply. (b): dropping the three state-gated candidates leaves 8, so the namespace is the capture
(the donor value is the capture's GSM) and the dropped captures are GSM5687824, GSM5687825 and
GSM5687826. Stated assumption, written into `reason`: one mouse per capture is assumed; the records do
not say which mice a capture holds or whether a capture pools mice. Daniel's Methods do not say either
(the one place they mention matching animals concerns the organ comparison, whose liver and lung
captures are out, and the declaration uses no clonal section). If
a capture pools mice, identical receptors from two mice there are one clone; if two captures share a
mouse, that sharing is not counted.

**Clone definition.** As `common.call_clones` implements it (spec §6.1): a contig qualifies when
`productive`, `high_confidence`, `full_length` and `is_cell` all parse true, `chain` is TRA or TRB and
`cdr3_nt` is ACGTN only; `clone_id` needs 1-2 distinct qualifying TRA and 1-2 TRB CDR3 nucleotide
sequences and is their sorted set; other cells are unlabelled. Clones are keyed by (namespace,
`clone_id`). `clone_id_trb` (1-2 TRB, the TRB set) is the sensitivity check. m = 2 primary, m = 3
sensitivity, as `[clonal]`.

**Root: `TEx_prog`.** Tpex is TCF-1+, Slamf6/Ly108+, PD-1 intermediate, TIM-3− (Im 2016: PD-1+
CXCR5+ Tim-3− TCF1+; Beltra 2020: TCF1+ PD-1 intermediate Ly108+; Utzschneider 2016, abstract:
Tcf1+ PD-1+). Daniel's Texprog is Tcf7+, Slamf6+, Id3+, sorted PD-1+SLAMF6+, and Daniel's introduction
defines Texprog as TCF1+CXCR5+PD-1+ citing Im and Utzschneider. Daniel does not state the PD-1 level
or TIM-3 for the cluster, so the match rests on TCF-1 and Slamf6; it is the only label carrying them
in exhausted cells. Corroboration from Task 3's counts: the SLAMF6+ sort's largest label is `TEx_prog`
(1,273 of 1,915). `TEx_prog` has 2,714 cells in the kept captures.

**Lineages.** End states that a paper shows arising from Tpex, and not converting into another
declared state, by a method other than clone sharing:
1. **Terminal exhausted (Texterm).** Hudson (CD101+Tim3+: arises from stem-like cells via the transitory
   cells, does not convert, transfer) and Beltra (Texterm: arises via Texint, keeps its state,
   transfer). Im 2016 shows CXCR5+Tim-3− giving rise to CXCR5−Tim-3+, a gate that Beltra shows also
   holds Texint cells (its Fig. S3B), so Im is not used for the end state. Zander does not qualify it: in
   Zander, CX3CR1−Ly108− donors partly acquire CX3CR1. The end maps to one label, `TEx`. The
   intermediate that Hudson and Beltra place on it (the transitory CD101−Tim3+ cells; Texint) is the
   CX3CR1+ subset and maps to **two labels**, `TEx_int` and `TEx_KLR` (section 2). Beltra's Texprog1
   and Texprog2 both map to the root.
2. **CX3CR1+ effector-like.** Zander: Ly108+ donors give rise to it, and its donors mostly keep their
   phenotype and become neither of the other subsets (transfer). Zander does not show it converting,
   so it is an end state by the rule, although Hudson and Beltra show the same subset converting to
   terminal cells in infection-matched hosts. Its end maps to **two labels**, `TEx_int` and `TEx_KLR`.

No other label is placed by a rule, so these are on no lineage: `Naive`, `Effector`,
`Effector_Memory`, `Memory`, `Ifn_Response`, `TEx_earlyEffector` (Daniel's own transfer of Texeeff
is not one of the three papers, and Texeeff precedes the root), and `TEx_lung` (a lung variant the
papers do not describe; the lung is out). TexKLR as a state of its own appears only in Daniel; its
end-state and bifurcation evidence (RNA velocity, Daniel's own transfer of Texint, and clones) is
not from Hudson, Zander or Beltra, so it places nothing. No link rests on an unreadable source.

**Registry form.** `RealDataset.lineages` holds each lineage as one label per stage, in order. A
stage that maps to two labels cannot be written that way without inventing an order between
`TEx_int` and `TEx_KLR`, or between `TEx_KLR` and `TEx`, which Daniel's own non-clonal transfer
(Texint gives rise to both TexKLR and Texterm) contradicts. Both lineages above have such a stage,
so `lineages = {}`. They are recorded here and in the registry's `notes`. A consequence: with no
lineage the order truth is all NaN, so no stage order is claimed that the literature does not give.

**Gated: none.** Gated means the root, intermediates and end each map to one label with more than 50
cells in the kept captures. The counts are no obstacle (kept-capture totals, from note section 3:
`TEx_prog` 2,714, `TEx_int` 7,495, `TEx_KLR` 5,606, `TEx` 8,403), but in both lineages one stage maps to
two labels, so `gated = []`.

**Two gated lineages: not assessable.** With 0 gated lineages there is no clonal comparison between
lineages, so `assessable = False`, reason "fewer than two gated lineages, so no clonal comparison
between lineages". The namespace rule (b) did not force not-assessable, so only this reason applies;
the namespace assumption is appended to it as the stated assumption. This outcome does not depend on
the two-label judgement above: mapping Zander's CX3CR1+ subset to `TEx_KLR` alone would give one
gated lineage, still fewer than two.

**Fallback.** Does not apply: `LABEL_COLUMN` is `"ident"`.

## 4. The declared lists

- `captures` (8): GSM5687827, GSM5687828, GSM5687829, GSM5687830, GSM5687831, GSM5687832,
  GSM5687833, GSM5687834. 33,188 cells by Task 3's counts.
- `namespace`: each capture maps to its own GSM.
- Excluded: acute GSM5687820-23 (arm); D5 GSM6381666-69 (separate matrix); liver and lung
  GSM5687816-19 (tissue); state-gated GSM5687824-26 (namespace rule b).
- `label_key = "state"`, `root_label = "TEx_prog"`, `lineages = {}`, `gated = []`.
- `clonal = RealClonal("clone_id", "donor", "clone_id_trb", "raw_clonotype_id", "capture", False,
  <reason>, "2026-11-20")`, with reason "fewer than two gated lineages, so no clonal comparison between
  lineages; namespace = capture (rule b, state-gated GSM5687824, GSM5687825 and GSM5687826 dropped):
  one mouse per capture is assumed; the records do not say which mice a capture holds or whether a
  capture pools mice".
- `report_time_key = "day"`; `time_key` stays unset (day is reported, never the order truth).
- Gates: `pt_spearman_lineage >= 0.70`, `branch_accuracy >= 0.80`, `lineage_recall >= 1.0`; the two
  clonal gates are left out (not assessable). `[real.daniel2022.clonal]` and
  `[real.daniel2022.reported]` as the plan gives them.

## 5. Unverified, and for the controller

- **Utzschneider 2016:** only the abstract was read; the root does not depend on it.
- **Mice per capture and pooling:** not stated by the records or by Daniel's Methods; assumed, as
  above.
- **Label names to clusters:** the name match is not stated as a table in the paper; it is read from
  the names and corroborated by Task 3's per-capture label counts.
- **Gap in the plan (not a rule conflict):** with `gated = []`, the three base gates score NaN, and
  `check_gates` skips NaN by default, so a Daniel row passes them without measuring anything. The
  "Two gated lineages" rule covers this for the clonal gates only. The gates are written as the plan
  gives them; Task 7/8 (or a controller ruling) should show "no gated lineage" beside every Daniel
  row so that a vacuous pass is not read as a pass.
