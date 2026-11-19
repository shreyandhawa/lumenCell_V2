# Daniel 2022 (GSE188666): sample records and metadata

Phase 3 Task 3, 2026-11-19. Written from the 46 GSM records of GSE188666 (GEO text form,
`acc.cgi?acc=<GSM>&targ=self&form=text&view=full`), their 46 linked BioSample records (NCBI
E-utilities `efetch db=biosample`, text), the GSE188666 and GSE188670 series records, and the two
metadata files, read by `benchmarks/inspection/daniel2022_metadata.py`. No matrix and no contig file
was opened. No page was behind a CAPTCHA. This note records what was found; the subset and the
namespace are Task 6's to declare.

Script runs (output kept under `$TMPDIR`, not committed):
- `uv run python -I benchmarks/inspection/daniel2022_metadata.py --dir <datasets dir>` → `daniel-metadata.txt`
- the same with `--labels ident --samples orig.ident` → `daniel-metadata-xtab.txt`
- the same with `--labels ident --samples lcmv`, and with `--labels tissue --samples lcmv`

## 1. Samples

46 GSMs: 23 gene-expression (GEX) libraries and 23 TCR libraries, one of each per capture.
**Every capture has a TCR file.** The spec's and the dataset note's "24 TCR files" was a miscount:
GSM5687816-GSM5687834 is 19 GSMs, not 20, and 19 + 4 = 23, which is what `GSE188666_RAW.tar` holds.

Fields, per column:
- **GSM**: `!Sample_geo_accession`.
- **Type**: `!Sample_description`, "5' gene expression library (10X Genomics)" (GEX) or "TCR library (10X Genomics)" (VDJ). `!Sample_library_strategy` is `RNA-Seq` for all 46, so it does not tell them apart.
- **Other GSM**: the GSM of the other type whose `!Sample_title` is the same apart from the `_gex`/`_tcr` ending. No `!Sample_relation` links them (each GSM's only relations are its own BioSample and SRA experiment).
- **Title**: `!Sample_title`.
- **Arm, day, tissue**: from `!Sample_source_name_ch1`, as written ("LCMV-c13" or "LCMV-Armstrong"; "<n> dpi"; "from <tissue>").
- **Sort gate**: `!Sample_characteristics_ch1`, `cell type: ...`, as written.
- **Mice**: no field of the GSM or BioSample record states the mice in a capture: `not stated` for all 46.
- **Shared mouse id**: `none`. The only per-sample identifiers are the GSM, the BioSample (46 distinct; their only attributes are `/source name` and `/cell type`, the same text as the GSM) and the SRA experiment.
- **Seq**: `!Sample_instrument_model` (N = Illumina NovaSeq 6000, H = Illumina HiSeq 4000), given because it bears on the pairing flags below.

| GSM | Type | Other GSM | Title | Arm | Day | Tissue | Sort gate | Mice | Shared mouse id | Seq |
|---|---|---|---|---|---|---|---|---|---|---|
| GSM5687797 | GEX | GSM5687816 | Liver_Chronic_D21_CD8_TetNeg_gex | LCMV-c13 | 21 dpi | liver | CD8+ gp33- T cells | not stated | none | N |
| GSM5687798 | GEX | GSM5687817 | Liver_Chronic_D21_CD8_TetPos_gex | LCMV-c13 | 21 dpi | liver | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687799 | GEX | GSM5687818 | Lung_Chronic_D21_CD8_TetNeg_gex | LCMV-c13 | 21 dpi | lung | CD8+ gp33- T cells | not stated | none | N |
| GSM5687800 | GEX | GSM5687819 | Lung_Chronic_D21_CD8_TetPos_gex | LCMV-c13 | 21 dpi | lung | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687801 | GEX | GSM5687820 | Spleen_Acute_D21_CD8_TetNeg_gex | LCMV-Armstrong | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687802 | GEX | GSM5687821 | Spleen_Acute_D21_CD8_TetPos_gex | LCMV-Armstrong | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687803 | GEX | GSM5687822 | Spleen_Acute_D8_CD8_TetNeg_gex | LCMV-Armstrong | 8 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687804 | GEX | GSM5687823 | Spleen_Acute_D8_CD8_TetPos_gex | LCMV-Armstrong | 8 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687805 | GEX | GSM5687825 | Spleen_Chronic_D21_CD8_CX3CR1Neg_gex | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1- T cells | not stated | none | H |
| GSM5687806 | GEX | GSM5687824 | Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos_gex | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1- SLAMF6+ T cells | not stated | none | H |
| GSM5687807 | GEX | GSM5687826 | Spleen_Chronic_D21_CD8_CX3CR1Pos_gex | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1+ T cells | not stated | none | H |
| GSM5687808 | GEX | GSM5687827 | Spleen_Chronic_D21_CD8_PD1Pos_gex | LCMV-c13 | 21 dpi | spleen | CD8+ PD1+ T cells | not stated | none | H |
| GSM5687809 | GEX | GSM5687828 | Spleen_Chronic_D21_CD8_TetMid_gex | LCMV-c13 | 21 dpi | spleen | CD8+ gp33 mid T cells | not stated | none | N |
| GSM5687810 | GEX | GSM5687830 | Spleen_Chronic_D21_CD8_TetNeg_gex | LCMV-c13 | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687811 | GEX | GSM5687829 | Spleen_Chronic_D21_CD8_TetNeg_rep1_gex | LCMV-c13 | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | H |
| GSM5687812 | GEX | GSM5687832 | Spleen_Chronic_D21_CD8_TetPos_gex | LCMV-c13 | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687813 | GEX | GSM5687831 | Spleen_Chronic_D21_CD8_TetPos_rep1_gex | LCMV-c13 | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | H |
| GSM5687814 | GEX | GSM5687833 | Spleen_Chronic_D8_CD8_TetNeg_gex | LCMV-c13 | 8 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687815 | GEX | GSM5687834 | Spleen_Chronic_D8_CD8_TetPos_gex | LCMV-c13 | 8 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687816 | VDJ | GSM5687797 | Liver_Chronic_D21_CD8_TetNeg_tcr | LCMV-c13 | 21 dpi | liver | CD8+ gp33- T cells | not stated | none | N |
| GSM5687817 | VDJ | GSM5687798 | Liver_Chronic_D21_CD8_TetPos_tcr | LCMV-c13 | 21 dpi | liver | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687818 | VDJ | GSM5687799 | Lung_Chronic_D21_CD8_TetNeg_tcr | LCMV-c13 | 21 dpi | lung | CD8+ gp33- T cells | not stated | none | N |
| GSM5687819 | VDJ | GSM5687800 | Lung_Chronic_D21_CD8_TetPos_tcr | LCMV-c13 | 21 dpi | lung | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687820 | VDJ | GSM5687801 | Spleen_Acute_D21_CD8_TetNeg_tcr | LCMV-Armstrong | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687821 | VDJ | GSM5687802 | Spleen_Acute_D21_CD8_TetPos_tcr | LCMV-Armstrong | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687822 | VDJ | GSM5687803 | Spleen_Acute_D8_CD8_TetNeg_tcr | LCMV-Armstrong | 8 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687823 | VDJ | GSM5687804 | Spleen_Acute_D8_CD8_TetPos_tcr | LCMV-Armstrong | 8 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687824 | VDJ | GSM5687806 | Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1- SLAMF6+ T cells | not stated | none | H |
| GSM5687825 | VDJ | GSM5687805 | Spleen_Chronic_D21_CD8_CX3CR1Neg_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1- T cells | not stated | none | H |
| GSM5687826 | VDJ | GSM5687807 | Spleen_Chronic_D21_CD8_CX3CR1Pos_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ CX3CR1+ T cells | not stated | none | H |
| GSM5687827 | VDJ | GSM5687808 | Spleen_Chronic_D21_CD8_PD1Pos_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ PD1+ T cells | not stated | none | H |
| GSM5687828 | VDJ | GSM5687809 | Spleen_Chronic_D21_CD8_TetMid_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ gp33 mid T cells | not stated | none | N |
| GSM5687829 | VDJ | GSM5687811 | Spleen_Chronic_D21_CD8_TetNeg_rep1_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687830 | VDJ | GSM5687810 | Spleen_Chronic_D21_CD8_TetNeg_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ gp33- T cells | not stated | none | H |
| GSM5687831 | VDJ | GSM5687813 | Spleen_Chronic_D21_CD8_TetPos_rep1_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM5687832 | VDJ | GSM5687812 | Spleen_Chronic_D21_CD8_TetPos_tcr | LCMV-c13 | 21 dpi | spleen | CD8+ gp33+ T cells | not stated | none | H |
| GSM5687833 | VDJ | GSM5687814 | Spleen_Chronic_D8_CD8_TetNeg_tcr | LCMV-c13 | 8 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM5687834 | VDJ | GSM5687815 | Spleen_Chronic_D8_CD8_TetPos_tcr | LCMV-c13 | 8 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM6381662 | GEX | GSM6381666 | Spleen_Chronic_D5_CD8_TetNeg_rep1_gex | LCMV-c13 | 5 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM6381663 | GEX | GSM6381667 | Spleen_Chronic_D5_CD8_TetPos_rep1_gex | LCMV-c13 | 5 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM6381664 | GEX | GSM6381668 | Spleen_Chronic_D5_CD8_TetNeg_rep2_gex | LCMV-c13 | 5 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM6381665 | GEX | GSM6381669 | Spleen_Chronic_D5_CD8_TetPos_rep2_gex | LCMV-c13 | 5 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM6381666 | VDJ | GSM6381662 | Spleen_Chronic_D5_CD8_TetNeg_rep1_tcr | LCMV-c13 | 5 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM6381667 | VDJ | GSM6381663 | Spleen_Chronic_D5_CD8_TetPos_rep1_tcr | LCMV-c13 | 5 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |
| GSM6381668 | VDJ | GSM6381664 | Spleen_Chronic_D5_CD8_TetNeg_rep2_tcr | LCMV-c13 | 5 dpi | spleen | CD8+ gp33- T cells | not stated | none | N |
| GSM6381669 | VDJ | GSM6381665 | Spleen_Chronic_D5_CD8_TetPos_rep2_tcr | LCMV-c13 | 5 dpi | spleen | CD8+ gp33+ T cells | not stated | none | N |

**Two record inconsistencies that Task 5a's join check must settle.** Neither can be settled from the records.
1. **D5 TCR titles against file names (`!Sample_supplementary_file_1`).** GSM6381667 is titled `..._D5_CD8_TetPos_rep1_tcr` and its source says gp33+, but its file is `GSM6381667_Spleen_Chronic_D5_CD8_TetNeg_rep2_tcr_filtered_contig_annotations.csv.gz`. GSM6381668 is titled `..._TetNeg_rep2_tcr` (gp33-), but its file is `GSM6381668_..._TetPos_rep1_tcr_...`. The other two D5 TCR GSMs agree with their files. The file names run in alphabetical order (TetNeg_rep1, TetNeg_rep2, TetPos_rep1, TetPos_rep2) while the titles follow the GEX order, which suggests the files were attached in the wrong order, but that is inference. So the capture of members `GSM6381667_*` and `GSM6381668_*` must come from the barcode shares, not from the member name or the GSM.
2. **D21 spleen tetramer replicates across sequencers.** Paired by title, GSM5687810 (TetNeg, NovaSeq) ↔ GSM5687830 (TetNeg, HiSeq), GSM5687811 (TetNeg_rep1, HiSeq) ↔ GSM5687829 (TetNeg_rep1, NovaSeq), and the same for TetPos (GSM5687812 N ↔ GSM5687832 H; GSM5687813 H ↔ GSM5687831 N). Every other pair has the same sequencer for its GEX and TCR libraries, and every other HiSeq library is in the D21 spleen CX3CR1/SLAMF6/PD-1 sorts. Paired by sequencer instead, "rep1" and the unnumbered replicate would swap. A GEX and a TCR library need not run on the same sequencer, so this is a flag, not a finding. The title pairing is used above; the join share of members GSM5687829-32 against the tags `..._TetNeg_gex_`, `..._TetNeg_rep1_gex_`, `..._TetPos_gex_` and `..._TetPos_rep1_gex_` decides it (those captures hold 8,914, 2,057, 8,194 and 2,031 cells in the main metadata).

## 2. Metadata

Both files are Seurat exports, tab-separated: the header has 11 fields and each row 12, the first being the cell id (pandas takes it as the index).

**`GSE188666_scRNA_LCMV_metadata.tsv.gz`**: 96,750 rows. First index values `Spleen_Chronic_D8_CD8_TetPos_gex_AAACCTGAGGCTAGAC-1`, `Spleen_Chronic_D8_CD8_TetPos_gex_AAACGGGCATTGTGCA-1`, `Spleen_Chronic_D8_CD8_TetPos_gex_AAACGGGTCACAGGCC-1`.

**`GSE188666_scRNA_D5_LCMV_metadata.txt.gz`**: 77,614 rows. First index values `Spleen_Chronic_D5_CD8_TetPos_rep2_gex_AAACCTGAGACTAAGT-1`, `Spleen_Chronic_D5_CD8_TetPos_rep2_gex_AAACCTGAGGACACCA-1`, `Spleen_Chronic_D5_CD8_TetPos_rep2_gex_AAACCTGAGGACATTA-1`.

Columns, the same in both files:

| Column | dtype | Withheld |
|---|---|---|
| orig.ident | str | - |
| nCount_RNA | int64 | - |
| nFeature_RNA | int64 | - |
| tcr_cdr3s_aa | str | withheld |
| tcr_cdr3s_nt | str | withheld |
| tissue | str | - |
| lcmv | str | - |
| sort | str | - |
| ident | str | - |
| UMAP_1 | float64 | - |
| UMAP_2 | float64 | - |

Other columns, as printed: `nCount_RNA` 505 / 4,500 / 128,904 (min / median / max) in the main file and 500 / 5,244 / 128,904 in the D5 file; `nFeature_RNA` 332 / 1,819 / 7,809 and 295 / 1,996 / 7,809. `tissue` is Liver, Lung or Spleen; `lcmv` is `<Arm>_D<day>` (`Acute_D8`, `Acute_D21`, `Chronic_D8`, `Chronic_D21`, and in the D5 file `Chronic_D5`); `sort` is TetNeg, TetPos, TetMid, CX3CR1Neg, CX3CR1Neg_Slamf6Pos, CX3CR1Pos or PD1Pos. `UMAP_1`/`UMAP_2` are each file's own embedding and are not used.

**Index tag vocabulary.** Every index value is `<tag prefix><16-nt core>-1`, the suffix always `-1`. The tag prefix is the GEX GSM title plus `_`, so each tag names one GEX GSM. **The sample column is `orig.ident`**: its value is the title without `_CD8` and `_gex`. Counts are cells per tag, which equal the cells per `orig.ident` value in every row.

| GEX GSM | Index tag | `orig.ident` | Cells, main file | Cells, D5 file |
|---|---|---|---|---|
| GSM5687797 | `Liver_Chronic_D21_CD8_TetNeg_gex_<core>-1` | Liver_Chronic_D21_TetNeg | 9,938 | - |
| GSM5687798 | `Liver_Chronic_D21_CD8_TetPos_gex_<core>-1` | Liver_Chronic_D21_TetPos | 3,847 | - |
| GSM5687799 | `Lung_Chronic_D21_CD8_TetNeg_gex_<core>-1` | Lung_Chronic_D21_TetNeg | 9,185 | - |
| GSM5687800 | `Lung_Chronic_D21_CD8_TetPos_gex_<core>-1` | Lung_Chronic_D21_TetPos | 5,670 | - |
| GSM5687801 | `Spleen_Acute_D21_CD8_TetNeg_gex_<core>-1` | Spleen_Acute_D21_TetNeg | 9,056 | - |
| GSM5687802 | `Spleen_Acute_D21_CD8_TetPos_gex_<core>-1` | Spleen_Acute_D21_TetPos | 7,654 | - |
| GSM5687803 | `Spleen_Acute_D8_CD8_TetNeg_gex_<core>-1` | Spleen_Acute_D8_TetNeg | 2,935 | - |
| GSM5687804 | `Spleen_Acute_D8_CD8_TetPos_gex_<core>-1` | Spleen_Acute_D8_TetPos | 4,621 | - |
| GSM5687805 | `Spleen_Chronic_D21_CD8_CX3CR1Neg_gex_<core>-1` | Spleen_Chronic_D21_CX3CR1Neg | 3,738 | 3,736 |
| GSM5687806 | `Spleen_Chronic_D21_CD8_CX3CR1Neg_Slamf6Pos_gex_<core>-1` | Spleen_Chronic_D21_CX3CR1Neg_Slamf6Pos | 1,915 | 1,915 |
| GSM5687807 | `Spleen_Chronic_D21_CD8_CX3CR1Pos_gex_<core>-1` | Spleen_Chronic_D21_CX3CR1Pos | 5,003 | 5,001 |
| GSM5687808 | `Spleen_Chronic_D21_CD8_PD1Pos_gex_<core>-1` | Spleen_Chronic_D21_PD1Pos | 1,989 | 1,987 |
| GSM5687809 | `Spleen_Chronic_D21_CD8_TetMid_gex_<core>-1` | Spleen_Chronic_D21_TetMid | 5,875 | 5,873 |
| GSM5687810 | `Spleen_Chronic_D21_CD8_TetNeg_gex_<core>-1` | Spleen_Chronic_D21_TetNeg | 8,914 | 8,904 |
| GSM5687811 | `Spleen_Chronic_D21_CD8_TetNeg_rep1_gex_<core>-1` | Spleen_Chronic_D21_TetNeg_rep1 | 2,057 | 2,057 |
| GSM5687812 | `Spleen_Chronic_D21_CD8_TetPos_gex_<core>-1` | Spleen_Chronic_D21_TetPos | 8,194 | 8,189 |
| GSM5687813 | `Spleen_Chronic_D21_CD8_TetPos_rep1_gex_<core>-1` | Spleen_Chronic_D21_TetPos_rep1 | 2,031 | 2,030 |
| GSM5687814 | `Spleen_Chronic_D8_CD8_TetNeg_gex_<core>-1` | Spleen_Chronic_D8_TetNeg | 2,526 | 2,526 |
| GSM5687815 | `Spleen_Chronic_D8_CD8_TetPos_gex_<core>-1` | Spleen_Chronic_D8_TetPos | 1,602 | 1,602 |
| GSM6381662 | `Spleen_Chronic_D5_CD8_TetNeg_rep1_gex_<core>-1` | Spleen_Chronic_D5_TetNeg_rep1 | - | 8,843 |
| GSM6381663 | `Spleen_Chronic_D5_CD8_TetPos_rep1_gex_<core>-1` | Spleen_Chronic_D5_TetPos_rep1 | - | 7,803 |
| GSM6381664 | `Spleen_Chronic_D5_CD8_TetNeg_rep2_gex_<core>-1` | Spleen_Chronic_D5_TetNeg_rep2 | - | 8,969 |
| GSM6381665 | `Spleen_Chronic_D5_CD8_TetPos_rep2_gex_<core>-1` | Spleen_Chronic_D5_TetPos_rep2 | - | 8,179 |

**The D5 file is not D5 only.** It holds the four D5 captures and the 11 chronic spleen D8 and D21 captures of the main file. Every one of its D8/D21 cell ids is also a main-file cell id (43,820 shared; 24 main-file chronic spleen cells are absent from it). This comparison used the index and `orig.ident`/`ident` only (`benchmarks/inspection/daniel2022_metadata.py --compare --labels ident --samples orig.ident`, which reads no other column); see HANDOFF "Phase 3 Task 3 decisions". Whether the D5 *matrix* holds the same D8/D21 barcodes is Task 5a's to check.

## 3. State labels

**The label column is `ident`.** The main file has 11 labels and the D5 file 7.

Main file, cells by arm and day (`lcmv`):

| Label | Acute_D8 | Acute_D21 | Chronic_D8 | Chronic_D21 |
|---|---|---|---|---|
| Effector | 6,338 | 39 | 294 | 707 |
| Effector_Memory | 76 | 4,640 | 7 | 558 |
| Ifn_Response | 78 | 15 | 102 | 1,018 |
| Memory | 105 | 7,506 | 25 | 1,003 |
| Naive | 57 | 4,096 | 323 | 5,204 |
| TEx | 1 | 11 | 9 | 24,912 |
| TEx_KLR | 392 | 341 | 430 | 11,683 |
| TEx_earlyEffector | 90 | 9 | 781 | 1,404 |
| TEx_int | 299 | 21 | 1,991 | 9,809 |
| TEx_lung | 8 | 2 | 59 | 7,928 |
| TEx_prog | 112 | 30 | 107 | 4,130 |

D5 file, cells by day (all chronic, all spleen):

| Label | Chronic_D5 | Chronic_D8 | Chronic_D21 |
|---|---|---|---|
| Naive | 4,319 | 315 | 3,963 |
| TEx | 1,207 | 88 | 13,087 |
| TEx_KLR | 2,844 | 551 | 8,202 |
| TEx_earlyEffector | 9,854 | 1,888 | 1,769 |
| TEx_earlyEffector2 | 2,184 | 318 | 528 |
| TEx_int | 10,418 | 833 | 7,201 |
| TEx_prog | 2,968 | 135 | 4,942 |

**The D5 metadata does not use the main vocabulary.** It shares six labels (Naive, TEx, TEx_KLR, TEx_earlyEffector, TEx_int, TEx_prog), adds `TEx_earlyEffector2`, and lacks Effector, Effector_Memory, Ifn_Response, Memory and TEx_lung. On the 43,820 cells present in both files the labels were assigned apart: for example, 2,316 cells labelled `TEx_int` in the main file are `TEx_earlyEffector` in D5, and 729 main-file `Naive` cells are `TEx_prog` in D5. Shared cells, main label (rows) by D5 label (columns):

| main \ D5 | Naive | TEx | TEx_KLR | TEx_earlyEffector | TEx_earlyEffector2 | TEx_int | TEx_prog |
|---|---|---|---|---|---|---|---|
| Effector | 0 | 0 | 386 | 202 | 0 | 253 | 20 |
| Effector_Memory | 0 | 0 | 438 | 0 | 0 | 16 | 1 |
| Ifn_Response | 0 | 692 | 10 | 10 | 0 | 39 | 8 |
| Memory | 70 | 28 | 44 | 6 | 1 | 66 | 294 |
| Naive | 4,191 | 27 | 18 | 1 | 13 | 9 | 729 |
| TEx | 2 | 11,036 | 0 | 19 | 0 | 246 | 14 |
| TEx_KLR | 0 | 5 | 7,626 | 34 | 0 | 755 | 3 |
| TEx_earlyEffector | 5 | 39 | 6 | 956 | 793 | 18 | 23 |
| TEx_int | 0 | 352 | 219 | 2,316 | 0 | 6,414 | 100 |
| TEx_lung | 1 | 964 | 2 | 89 | 1 | 201 | 8 |
| TEx_prog | 9 | 32 | 4 | 24 | 38 | 17 | 3,877 |

Main file, cells per capture (`orig.ident`, rows) and label (columns), for Task 6's "more than 50 cells in the kept captures":

| orig.ident | Effector | Effector_Memory | Ifn_Response | Memory | Naive | TEx | TEx_KLR | TEx_earlyEffector | TEx_int | TEx_lung | TEx_prog |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Liver_Chronic_D21_TetNeg | 7 | 7 | 104 | 43 | 122 | 9216 | 103 | 21 | 149 | 140 | 26 |
| Liver_Chronic_D21_TetPos | 0 | 0 | 32 | 1 | 0 | 3661 | 13 | 4 | 62 | 65 | 9 |
| Lung_Chronic_D21_TetNeg | 120 | 97 | 170 | 465 | 417 | 465 | 2658 | 165 | 1228 | 3265 | 135 |
| Lung_Chronic_D21_TetPos | 13 | 6 | 53 | 10 | 0 | 247 | 916 | 155 | 953 | 3251 | 66 |
| Spleen_Acute_D21_TetNeg | 31 | 2132 | 1 | 2615 | 4035 | 7 | 193 | 5 | 18 | 2 | 17 |
| Spleen_Acute_D21_TetPos | 8 | 2508 | 14 | 4891 | 61 | 4 | 148 | 4 | 3 | 0 | 13 |
| Spleen_Acute_D8_TetNeg | 2282 | 58 | 23 | 61 | 52 | 1 | 217 | 23 | 178 | 6 | 34 |
| Spleen_Acute_D8_TetPos | 4056 | 18 | 55 | 44 | 5 | 0 | 175 | 67 | 121 | 2 | 78 |
| Spleen_Chronic_D21_CX3CR1Neg | 12 | 2 | 92 | 44 | 10 | 2893 | 12 | 31 | 293 | 340 | 9 |
| Spleen_Chronic_D21_CX3CR1Neg_Slamf6Pos | 11 | 1 | 9 | 212 | 163 | 16 | 0 | 121 | 99 | 10 | 1273 |
| Spleen_Chronic_D21_CX3CR1Pos | 287 | 169 | 42 | 31 | 2 | 20 | 2805 | 65 | 1521 | 56 | 5 |
| Spleen_Chronic_D21_PD1Pos | 59 | 25 | 43 | 45 | 16 | 563 | 469 | 29 | 520 | 78 | 142 |
| Spleen_Chronic_D21_TetMid | 43 | 35 | 85 | 30 | 701 | 1471 | 1785 | 163 | 951 | 152 | 459 |
| Spleen_Chronic_D21_TetNeg | 56 | 97 | 138 | 60 | 2936 | 2006 | 1385 | 224 | 1311 | 225 | 476 |
| Spleen_Chronic_D21_TetNeg_rep1 | 40 | 83 | 11 | 28 | 812 | 457 | 171 | 71 | 199 | 86 | 99 |
| Spleen_Chronic_D21_TetPos | 26 | 4 | 208 | 12 | 19 | 2920 | 1169 | 296 | 2134 | 214 | 1192 |
| Spleen_Chronic_D21_TetPos_rep1 | 33 | 32 | 31 | 22 | 6 | 977 | 197 | 59 | 389 | 46 | 239 |
| Spleen_Chronic_D8_TetNeg | 228 | 7 | 69 | 22 | 321 | 6 | 414 | 240 | 1137 | 29 | 53 |
| Spleen_Chronic_D8_TetPos | 66 | 0 | 33 | 3 | 2 | 3 | 16 | 541 | 854 | 30 | 54 |

D5 file, the same by capture:

| orig.ident | Naive | TEx | TEx_KLR | TEx_earlyEffector | TEx_earlyEffector2 | TEx_int | TEx_prog |
|---|---|---|---|---|---|---|---|
| Spleen_Chronic_D21_CX3CR1Neg | 9 | 3241 | 19 | 90 | 10 | 334 | 33 |
| Spleen_Chronic_D21_CX3CR1Neg_Slamf6Pos | 65 | 17 | 2 | 96 | 111 | 52 | 1572 |
| Spleen_Chronic_D21_CX3CR1Pos | 2 | 78 | 2793 | 398 | 6 | 1696 | 28 |
| Spleen_Chronic_D21_PD1Pos | 7 | 675 | 452 | 136 | 12 | 512 | 193 |
| Spleen_Chronic_D21_TetMid | 624 | 1693 | 1831 | 160 | 70 | 947 | 548 |
| Spleen_Chronic_D21_TetNeg | 2489 | 2465 | 1505 | 116 | 152 | 1261 | 916 |
| Spleen_Chronic_D21_TetNeg_rep1 | 749 | 536 | 270 | 114 | 39 | 161 | 188 |
| Spleen_Chronic_D21_TetPos | 14 | 3350 | 1108 | 481 | 107 | 1922 | 1207 |
| Spleen_Chronic_D21_TetPos_rep1 | 4 | 1032 | 222 | 178 | 21 | 316 | 257 |
| Spleen_Chronic_D5_TetNeg_rep1 | 2592 | 547 | 1138 | 1300 | 321 | 2313 | 632 |
| Spleen_Chronic_D5_TetNeg_rep2 | 1675 | 399 | 1255 | 1584 | 463 | 2554 | 1039 |
| Spleen_Chronic_D5_TetPos_rep1 | 22 | 137 | 215 | 3290 | 712 | 2690 | 737 |
| Spleen_Chronic_D5_TetPos_rep2 | 30 | 124 | 236 | 3680 | 688 | 2861 | 560 |
| Spleen_Chronic_D8_TetNeg | 313 | 51 | 523 | 912 | 74 | 576 | 77 |
| Spleen_Chronic_D8_TetPos | 2 | 37 | 28 | 976 | 244 | 257 | 58 |

Tissue by arm and day (main file): liver and lung appear only at Chronic_D21 (13,785 and 14,855 cells); spleen at every arm and day (Acute_D8 7,556, Acute_D21 16,710, Chronic_D8 4,128, Chronic_D21 39,716). `TEx_lung` is a tissue-named label: 6,516 of its 7,997 cells are in the two lung captures.

## 4. Sort gates

From `!Sample_characteristics_ch1` (`cell type:`). Each gate is listed with the GEX and TCR GSMs that use it; the metadata `sort` value is given in brackets.

| Gate (as written) | Markers | `sort` | GEX GSMs | TCR GSMs |
|---|---|---|---|---|
| CD8+ gp33- T cells | CD8; gp33 negative | TetNeg | GSM5687797, GSM5687799, GSM5687801, GSM5687803, GSM5687810, GSM5687811, GSM5687814, GSM6381662, GSM6381664 | GSM5687816, GSM5687818, GSM5687820, GSM5687822, GSM5687829, GSM5687830, GSM5687833, GSM6381666, GSM6381668 |
| CD8+ gp33+ T cells | CD8; gp33 positive | TetPos | GSM5687798, GSM5687800, GSM5687802, GSM5687804, GSM5687812, GSM5687813, GSM5687815, GSM6381663, GSM6381665 | GSM5687817, GSM5687819, GSM5687821, GSM5687823, GSM5687831, GSM5687832, GSM5687834, GSM6381667, GSM6381669 |
| CD8+ gp33 mid T cells | CD8; gp33 intermediate | TetMid | GSM5687809 | GSM5687828 |
| CD8+ CX3CR1- T cells | CD8; CX3CR1 negative | CX3CR1Neg | GSM5687805 | GSM5687825 |
| CD8+ CX3CR1- SLAMF6+ T cells | CD8; CX3CR1 negative; SLAMF6 positive | CX3CR1Neg_Slamf6Pos | GSM5687806 | GSM5687824 |
| CD8+ CX3CR1+ T cells | CD8; CX3CR1 positive | CX3CR1Pos | GSM5687807 | GSM5687826 |
| CD8+ PD1+ T cells | CD8; PD-1 positive | PD1Pos | GSM5687808 | GSM5687827 |

"gp33" is the LCMV GP33-41 epitope; that the gate is a tetramer stain is read from the titles (`TetPos`, `TetNeg`, `TetMid`) and `!Sample_series_overall_design` ("sorted to obtain either gp-33 positive or negative CD8+ T cells"); the tetramer's MHC allele is not stated. The records do not say whether the CX3CR1, SLAMF6 and PD-1 sorts were also tetramer-gated. The two GSM6381667/GSM6381668 entries follow the titles; see the flag in section 1.

## 5. Per-mouse identifiability

- **Is a mouse recorded for each capture? No.** No GSM field (title, source, characteristics, description, relations) and no BioSample attribute names a mouse, an animal number or a pool. The series' overall design does not either.
- **Does any mouse span several captures? Not stated.** No identifier is shared between GSMs except the series. Paired GEX and TCR libraries share a title but have separate BioSamples, so even that pairing is by title only. `!Sample_series_summary` says the paper compares "shared T cell clones that traffic to multiple organs", which implies that some liver, lung and spleen captures came from the same mice, but the records do not say which; this is inference, unverified, and rests on the authors' clonal analysis.
- **Does any capture pool mice? Not stated.** No record mentions pooling, hashing or a mouse count. Captures named `rep1`/`rep2` (D21 spleen TetNeg and TetPos; all four D5) are replicates of a sort, but the records do not say whether a replicate is another mouse.

## 6. Processing

**Not stated.** Each GSM says only "Processed with cellranger pipeline, version 3.1.0" and genome mm10. Nothing says the main and D5 matrices came from one Cell Ranger run, one `aggr`, or one integration. Evidence that they were handled apart:
- The four D5 GEX and four D5 TCR GSMs were submitted on Jul 25 2022; the other 38 on Nov 12 2021.
- The 38 earlier records name the main files (`Supplementary_files_format_and_content: scRNA_LCMV.mtx, ...`); the D5 records name no processed file and write "Assembly: mm10" where the others write "Genome_build: mm10".
- The D5 metadata is a separate Seurat object: its own cells (it repeats 43,820 of the main file's chronic spleen D8/D21 cells and drops 24), its own embedding and its own labels (section 3).

## 7. Withheld

Withheld by the freeze rule, in both metadata files: `tcr_cdr3s_aa` (str) and `tcr_cdr3s_nt` (str). They were listed by name and dtype only; no value was printed, summarised or used.

No clone × state, clone × capture or clone × lineage count was computed.

## 8. Still unverified

- **The paper's Methods.** Mice per capture, pooling, hashing and how the two Seurat objects were made are not checked: PMC11225711 is behind a CAPTCHA (Task 2, not bypassed) and Europe PMC serves no full text for it (it is not open access; `fullTextXML` returned HTTP 500 on 2026-11-19).
- **Which TCR member belongs to which capture,** for the two flags in section 1 (GSM6381667/GSM6381668, and GSM5687829-GSM5687832). Task 5a's join check decides.
- **What the D5 matrix holds.** Whether its barcodes are the D5 file's 77,614 cells, including the repeated D8/D21 cells, is Task 5a's to check.
- **The gate markers' reagents:** the gp33 tetramer's allele, and whether the CX3CR1, SLAMF6 and PD-1 sorts were tetramer-gated.
- **What `ident` labels mean.** Their marker definitions are in the paper, which was not read here; Task 6 maps labels to subsets from the paper's non-clonal marker descriptions.
