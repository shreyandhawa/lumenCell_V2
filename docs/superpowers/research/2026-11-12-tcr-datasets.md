# TCR + expression datasets for lumenCell: verification notes (2026-11-12)

Method: GEO text records (acc.cgi ... targ=self&form=text&view=full), NCBI FTP directory listings and filelist.txt (names and sizes only), PMC and Europe PMC pages. No data file was downloaded. The normal GEO HTML pages and the NCBI disclaimer page returned a reCAPTCHA wall; I did not bypass it. Anything not read from a primary source is marked UNVERIFIED. Nature.com article pages redirect to an auth gate, so paper text came from PMC.

## Verdicts

| Role | Candidate | Verdict |
|---|---|---|
| Benchmark | Daniel 2022, GSE188670 | Usable with caveats |
| Case study | Sheih 2020, GSE125881 | Usable with caveats, weaker than the alternative. TCR is in a non-standard "unpaired clonotypes" format, and cell-level joinability is UNVERIFIED |
| Case study alternative | Cao 2025, GSE207192 (axi-cel, DLBCL) | Usable with caveats. Better TCR format, but a non-commercial/no-derivatives article licence and only 7 patients |

## 1. Benchmark: Daniel et al. 2022 (GSE188670)

### 1.1 Accession matches the paper: YES
- Paper: Daniel B, Yost KE, ... Satpathy AT. "Divergent clonal differentiation trajectories of T cell exhaustion". Nature Immunology 23(11):1614-1627, 2022 (epub 26 Oct 2022). DOI 10.1038/s41590-022-01337-5. PMID 36289450, PMCID PMC11225711.
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC11225711/
  - https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:36289450&resultType=core&format=json
- GEO GSE188670 has the same title. It is a SuperSeries, not a data series. Public 25 Jul 2022, last updated 9 Feb 2023, Mus musculus, 57 GSMs, contact A. Satpathy (Stanford), BioProject PRJNA779919, PMID 36289450.
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188670 (read as text form)
- The paper's Data availability states "GEO accession: GSE188670". Code: https://github.com/katieyost/LCMV-code-2022 (repo contents not inspected, UNVERIFIED).
- SubSeries (the SuperSeries has no separate expression files of its own):
  - GSE188666 = scRNA-seq + scTCR-seq (46 GSMs: GSM5687797-GSM5687834 and GSM6381662-GSM6381669; PRJNA779914; SRP345737). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188666
  - GSE188669 = scATAC-seq (11 GSMs; PRJNA779923; SRP345739). https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE188669

### 1.2 Design
- Organism: mouse. LCMV: Clone 13 (chronic) or Armstrong (acute). CD8+ T cells sorted at day 5, 8 or 21 post infection (GSE188666 overall design).
- Tissues (from file names): spleen, liver, lung. Liver and lung are chronic D21 only. Acute is spleen D8 and D21.
- Sorting (file names plus paper): gp33 tetramer Pos / Neg / Mid, PD1+, CX3CR1+/-, Slamf6+. Cells are FACS-sorted on a BD FACSAria III.
- Chemistry (paper Methods, PMC): "10X 5' Single Cell Immune Profiling Solution Kit (v1.1 Chemistry)" for expression, "Chromium Single Cell V(D)J Enrichment Kit for mouse T cells (v1.1)" for TCR. About 10,000 cells per capture. cellranger count/vdj v3.1.0. No 3' or v2 chemistry.
- Cell numbers (paper Results as relayed from PMC): 96,750 scRNA profiles passed QC; TCR detected in 88,696 (91.7%); 62,731 scATAC profiles. Per-sample counts: UNVERIFIED (in the metadata file, not read).
- Mice: male C57BL/6J (CD45.2) and CD45.1, 8-12 wk. Number of mice per sequencing sample is NOT stated in Methods. No hashing or pooled-mouse deconvolution is mentioned. UNVERIFIED whether a capture pools several mice. Consequence: "donor" cannot be resolved below sample/capture. Treat each GSM (capture) as the clone-label namespace, which is the safe choice. Same CDR3 in two different mice or captures is not evidence of a shared ancestor.
- Sample count: GSE188666 lists 46 GSMs. 24 of them have a per-sample TCR file (20 in GSM5687816-34 plus 4 day-5 GSMs GSM6381666-69). The remaining GSMs presumably are the matching expression libraries; the mapping of GSM to the pooled matrices is UNVERIFIED.

### 1.3 Files and sizes
GSE188670 suppl (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188670/suppl/, filelist.txt): GSE188670_RAW.tar = 14,670,489,600 bytes (14 GB), which is the 20+4 TCR csv.gz plus 11 ATAC fragment files. Do not download this one.

GSE188666 suppl (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/):
| File | Size |
|---|---|
| GSE188666_scRNA_LCMV.mtx.gz | 574 MB |
| GSE188666_scRNA_LCMV_barcodes.tsv.gz | 476 KB |
| GSE188666_scRNA_LCMV_features.tsv.gz | 111 KB |
| GSE188666_scRNA_LCMV_metadata.tsv.gz | 4.6 MB |
| GSE188666_scRNA_D5_LCMV.mtx.gz | 519 MB |
| GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz | 381 KB |
| GSE188666_scRNA_D5_LCMV_features.tsv.gz | 111 KB |
| GSE188666_scRNA_D5_LCMV_metadata.txt.gz | 3.8 MB |
| GSE188666_RAW.tar | 12,113,920 bytes (about 12 MB; 24 TCR csv.gz) |

- Expression: 10x-style MatrixMarket + barcodes + features, two merged objects (main set, and day 5 set added later). Raw vs normalised counts: UNVERIFIED (not stated on GEO page; the files are .mtx, the paper used Seurat, so the matrix may be counts; check the first lines of the mtx header or the metadata before relying on it).
- No .h5ad and no .rds on GEO. Only mtx/tsv.
- Per-sample TCR: 24 files `GSM..._tcr_filtered_contig_annotations.csv.gz` (10x VDJ filtered_contig_annotations, the standard format, sizes 114 KB to 1.19 MB each). Names carry tissue, infection, day, sort gate. No joined obs column on GEO; the metadata.tsv column list is UNVERIFIED.
- Joining: the TCR files carry 10x barcodes per sample (standard, UNVERIFIED by inspection). The merged barcodes file probably carries a sample suffix or prefix; how it maps to the GSM names is UNVERIFIED. This is the first thing to inspect (read just the head of the barcodes and metadata files) before building on it.
- Download totals:
  - Minimum useful set (both matrices, barcodes, features, metadata, plus GSE188666_RAW.tar): about 574+519+0.5+0.4+0.2+4.6+3.8+12 = about 1.11 GB.
  - Whole SuperSeries RAW.tar alone: 14.7 GB (not needed).
  - The 24 TCR files are also available as individual files via the FTP suppl listing (sum of listed sizes about 12 MB).

### 1.4 Licence and terms
- Open public GEO series; no controlled access (mouse, not human, and in standard GEO). 
- GEO/NCBI policy page could not be read (reCAPTCHA): UNVERIFIED verbatim. General NLM policy read at https://www.nlm.nih.gov/web_policies.html says non-government content may be copyrighted and reuse is the user's responsibility; GEO's own statement that submitted data carry no use restriction was not confirmed from the GEO page. 
- Paper licence: PMC page shows only generic copyright and disclaimer links; no CC licence stated on the page I read. UNVERIFIED. Treat as: cite the paper and the accession.
- Code repo licence: UNVERIFIED.

### 1.5 Truth definition options
Independent of trajectory inference:
1. Day post infection (5/8/21) from the sample name: experimental, not inferred. Day 5 is early and day 21 late. Use it only within one infection arm (chronic vs acute) and tissue, since sorting differs by arm.
2. Acute vs chronic arm and tissue (liver, lung vs spleen): experimental labels (liver is reported to drive stronger exhaustion; treat as a biological covariate, not an ordering).
3. Sorting gates (tetramer Pos/Neg/Mid, PD1, CX3CR1, Slamf6): experimental but used with other info; they are biased to cell states, so they are partial truth (e.g. Slamf6+CX3CR1- should be progenitor-like, CX3CR1+ later), based on known biology, not on this dataset's trajectory.
4. Paper's cluster labels (Tnaive, Teff, Tem, Tmem, Texprog, Texint, Texterm, Texeeff, TexKLR, Texlung, TexISG): derived from expression clustering and marker genes. Not trajectory-derived, but not independent of expression, so they are circular for expression-based TI. Good for naming end states and the root; weak as independent evidence.
Not independent: the paper's pseudotime (dynamo UMAP potential, velocyto, ArchR addTrajectory). Do not use as truth.
TCR-based truth: the paper's claim of terminal-biased vs KLR-biased vs divergent clones is itself the finding being tested; it relies on clonotypes. Using it as truth would be circular with lumenCell's TCR evidence.
Convention suggestion for lumenCell: time point as ordering truth, cluster labels for end-state identity, TCR left as the independent evidence.

## 2. Case study: Sheih et al. 2020 (GSE125881)

### 2.1 Accession matches the paper: YES, title differs
- Paper: Sheih A, Voillet V, Hanafi L-A, ... Gottardo R, Linsley PS, Turtle CJ. "Clonal kinetics and single-cell transcriptional profiling of CAR-T cells in patients undergoing CD19 CAR-T immunotherapy". Nature Communications 11:219, 10 Jan 2020. DOI 10.1038/s41467-019-13880-1. PMID 31924795 (from the GEO record), PMCID PMC6954177. Licence CC BY 4.0.
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC6954177/
- GEO GSE125881 title: "Clonal kinetics and single-cell transcriptional profiling of adoptively transferred CD19-specific CD8+ CAR-T cells in adults with B cell malignancies". Not identical to the paper's title but the same study (same PMID, authors, data availability says GSE125881). Submitted 30 Jan 2019, public 22 Nov 2019, updated 21 Jan 2020. Homo sapiens, GPL16791, BioProject PRJNA517816, SRA SRP182902.
  - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125881 (text form)
- The accession was "believed" correct; it is correct.

### 2.2 Design
- 4 patients: CLL-1, CLL-2, NHL-6, NHL-7 (from GSM file names). The paper says 2 NHL and 2 CLL. 16 time-point libraries, 4 per patient: IP (infusion product), early (d12/d21), late (d28/d29/d38) and very late (d83-d112). Exact days from file names: CLL-1 IP,d21,d38,d112; CLL-2 IP,d12,d29,d83; NHL-6 IP,d12,d29,d102; NHL-7 IP,d12,d28,d89.
- GEO lists 32 GSMs (GSM3583861-92): 16 gene-expression + 16 VDJ.
- CD8+ CAR-T cells only, sorted CD8+/CD4-/EGFRt+ on BD Aria 2. Product: TCM-derived (n=2) or bulk CD8+ (n=2) CAR-T with 4-1BB/CD3z CAR; trial NCT01865617 (paper text has an NCT03690193 link discrepancy; not important).
- Tissue: blood after infusion, and the infusion product.
- Cells: 62,167 CD8+ CAR-T cells after QC (2,340 removed). Per-sample counts: UNVERIFIED.
- Chemistry: Chromium Single Cell 5' Gene Expression and V(D)J enrichment, HiSeq 2500, Cell Ranger 2.1.0, STAR, GRCh38 (so v1 era chemistry, version UNVERIFIED), aggregated with cellranger aggr.

### 2.3 Files
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE125nnn/GSE125881/suppl/
| File | Size |
|---|---|
| GSE125881_raw.expMatrix.csv.gz | 109 MB |
| GSE125881_RAW.tar | 2.0 MB (16 clonotype files) |
| 16 x GSM3583877-92_*-unpaired-clonotypes.csv.gz | 6.9 KB to 297 KB each |
- Total about 111 MB. A single aggregated UMI count matrix (csv, "raw" per file name and GSM text). No .h5ad, no .rds, no 10x mtx. A csv of 62k cells x genes at 109 MB gz is plausible only if sparse-ish; its orientation (genes x cells) and barcode format: UNVERIFIED.
- Sample encoding: GSM text says "Sample ID in raw.expmatrix.csv: 1" for CLL-1 IP, i.e. the aggr sample index (the 10x barcode suffix -1..-16 presumably). UNVERIFIED by looking at the matrix.
- TCR: provided as "unpaired clonotype tables" (GSM text: "V(D)J data were provided as unpaired clonotype tables"). That is NOT the standard filtered_contig_annotations. Whether the rows are per cell (with barcode) or per clonotype: UNVERIFIED. If there is no barcode column, a per-cell join is impossible from GEO; raw reads are in SRA (SRP182902) and would need re-running cellranger vdj. GEO also notes that NHL-6-IP and NHL-7-IP shared a chromium index, so their VDJ barcodes were used to tell them apart, which suggests barcodes do exist in the VDJ tables (inference, UNVERIFIED).
- Code: https://github.com/ValentinVoillet/CAR-T (MIT; R_code.R, meta-data.csv; Seurat v2 code). It may show how TCR is linked to cells (not inspected beyond file list).

### 2.4 Licence and access
- Public GEO, human but not controlled-access: the paper's data availability lists GEO GSE125881 for scRNA-seq, SRA PRJNA589633 for integration sites, and Source Data for TCRB. No dbGaP or controlled-access statement was seen (I read the first 100,000 of about 117,000 characters). Whether patient-level SRA reads are open or restricted: UNVERIFIED (GEO lists SRP182902).
- Paper CC BY 4.0. GEO policy: UNVERIFIED (see 1.4).

### 2.5 Truth for a case study (no benchmark truth needed)
- Time point (IP -> early -> late -> very late) per patient is independent of TI. The four IP clusters (cytotoxic, proliferating, pre-apoptotic mitochondrial, activation) are expression-derived. The paper's clonal kinetics (IRF/DRF clones) come from NHL-6 and NHL-7 only.

### 2.6 Verdict
Usable with caveats: correct accession, open, small (about 111 MB), CC BY. But TCR is "unpaired clonotypes", not 10x VDJ contig annotations, and cell-level joinability is not verified. Confirm barcode column and aggr suffix before committing; otherwise fall back to the alternative.

## 3. Case-study alternative: Cao et al. 2025 (GSE207192)

- Paper: Cao G, Hu Y, Pan T, Tang E, Asby N, Althaus T, Wan J, Riedell PA, Bishop MR, Kline JP, Huang J. "Two-stage CD8+ CAR T-cell differentiation in patients with large B-cell lymphoma". Nat Commun 16:4205 (May 2025). DOI 10.1038/s41467-025-59298-w. PMID 40328775, PMCID PMC12055983. Licence CC BY-NC-ND 4.0.
  - https://pmc.ncbi.nlm.nih.gov/articles/PMC12055983/
  - https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:40328775&resultType=core&format=json
- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207192 (text form). Title "Two-Stage CAR T-Cell Differentiation in Patients with Large B-Cell Lymphoma". Human, 60 GSMs, public 3 Apr 2025, BioProject PRJNA854096. Note the GEO design text says 3 patients while the paper says 7 patients (P1-P7) and about 38k CD8 cells (~67k CD4+CD8 in Table 2). Conflict between GEO and the paper is UNRESOLVED; GEO sample counts must be reconciled with the patient identifiers in metadata.
- Design: axi-cel (CD28 costimulated CD19 CAR), DLBCL, complete responders. Time points: IP, Texp (d8-14), Tper1 (d21), Tper2 (d28; P3 lacks Tper2, P6 lacks IP in Table 2). Blood PBMC and IP; CD19-tetramer sorting. Chemistry: Chromium Next GEM 5' Kit v2 with TotalSeq-C CITE-seq (RNA + ADT) and Human TCR Amplification Kit; NovaSeq 6000/NextSeq 550; Cell Ranger 7.0.0 (expression), 6.0.0 (vdj).
- Files (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207192/suppl/), 18 files:
  - TCR: All_TCR_Aggregated_clonotypes.csv.gz 1.7 MB; ..._consensus_annotations.csv.gz 11 MB; ..._filtered_contig_annotations.csv.gz 16 MB (the standard 10x format, aggregated across samples).
  - Post-infusion: scRNA_Raw_Counts.mtx.gz 294 MB, ADT mtx 4.1 MB, features, cells_colnames 255 KB, metaData 319 KB.
  - Pre-infusion: scRNA_Raw_Counts.mtx.gz 245 MB, ADT mtx 2.3 MB, cells_colnames 104 KB, metaData 134 KB.
  - CD4 CAR-T: counts RNA 6.6 MB, ADT 105 KB, metadata 22 KB.
  - Total about 590 MB. Raw counts, mtx, no h5ad or rds on GEO (an "R object of combined scTCR" is mentioned in the GSM text but I see no .rds in the listing; UNVERIFIED).
- Joining: the GSM text says Cell Ranger aggr was run for TCR with donor and timepoint as metadata, so aggregated contig barcodes carry a numeric suffix. Mapping of that suffix to the cells_colnames files and metadata (which has patient/timepoint): UNVERIFIED. Clonotypes were defined by paired productive TCRa/TCRb CDR3 nt, and the authors re-grouped clones within a patient, so keep clone labels per patient as required.
- Access: public GEO; no controlled-access statement; the paper does not say access is controlled (read first 100,000 of 143,722 characters). Licence CC BY-NC-ND for the article; GEO terms UNVERIFIED. The ND clause concerns article text, not obviously the data; UNVERIFIED legal reading.
- Verdict: usable with caveats. Better TCR format than Sheih, with multi-timepoint design and a paper that is explicitly about clone-state dynamics. Sample counts conflict (GEO says 3 patients; paper says 7).

## 4. Other alternatives noted (not verified, search results only)
- GSE206230 (Kasmani/Zander, JEM 2023, chronic LCMV Clone 13 lineage tracing with paired scRNA + scTCR). Not opened. UNVERIFIED. Source: https://rupress.org/jem/article/220/1/e20220679/213624/Clonal-lineage-tracing-reveals-mechanisms-skewing
- LRT (PLOS Comp Biol) applies acute LCMV paired scRNA/scTCR; accession not found. https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1011300
- The benchmark candidate passed, so no further mouse alternative was verified.

## 5. Key risks to carry into design
1. Daniel: the donor (mouse) is not identifiable; use the GSM capture as the clone namespace. Confirm barcode-to-GSM mapping from the merged barcodes/metadata files first.
2. Daniel: day 5 lives in a separate matrix from days 8 and 21; the two were not necessarily integrated together.
3. Sheih: TCR table format and barcode join unverified; if there are no barcodes, the dataset fails the "TCR per cell" requirement on GEO and needs SRA re-processing.
4. Cluster labels are not independent of expression, and the paper's pseudotime is not independent of trajectory inference.
5. GEO/NCBI terms were not read from the primary page (reCAPTCHA). Re-check https://www.ncbi.nlm.nih.gov/home/about/policies/ and https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html in a browser.
