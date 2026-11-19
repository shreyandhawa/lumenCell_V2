# Sheih 2020 (GSE125881): do the TCR tables join to cells? (2026-11-19)

**Verdict: PASS.** Sheih 2020 stays the case study; the Cao 2025 branch is not taken and nothing was downloaded from GSE207192.

Pass rule (fixed in the plan, Task 9, before looking): the clonotype table carries per-cell barcodes (a column, or lists per clonotype), and at least 0.5 of its distinct barcodes are among the matrix's cell IDs under one stated mapping.

## Result

- Table checked: `GSM3583882_CLL-2_d12-unpaired-clonotypes.csv.gz`, the largest post-infusion table (305,886 bytes in the tar, the largest of all 16 members; the largest infusion-product table is 296,504). Owner OK: Shreyan, in chat on 2026-11-19, "Yes, approve all" (relayed by the parent session).
- Barcode column: `barcodes` (chosen as the column with the most distinct barcodes among those where more than half of the non-null cells hold only barcodes; the script exits 1 if there is none). Each row is one clonotype x chain; the column holds `;`-separated lists of 10x barcodes (`[ACGT]{16}-1`; up to 23 per row; 1,579 of 8,391 rows hold more than one). No `,` separators.
- **Mapping rule:** a table barcode `<16-nt core>-1` maps to the matrix cell ID `<16-nt core>-<aggr sample ID>`, where the aggr sample ID is the GSM's "id in raw.expmatrix.csv".
- **Mapping details:** only `-1` barcodes join; a table barcode with any other ending, bare 16-nt ones included, does not join, stays in the denominator and is counted in the output (here 0 of 5,873). The share's denominator is the distinct table barcodes.
- **Negative control:** the same table mapped with the wrong aggr suffix `-10` against the same 6,271 `-9` cell IDs: 0 of 5,873, share 0.00.
- Matrix cell IDs with suffix `-9` (aggr ID of GSM3583882): 6,271. Distinct table barcodes: 5,873. **Found among the matrix cell IDs: 5,817 of 5,873, share 0.99** (pass rule: at least 0.50).
- Informational, same rule on the other 15 tables (every one has at least 0.97 of its distinct barcodes in the matrix; the lowest are the infusion-product tables of NHL-6 and NHL-7, which share a chromium index): aggr 1: 5272/5346, 2: 6128/6140, 3: 3257/3265, 4: 718/718, 8: 1659/1678, 10: 4021/4028, 11: 548/548, 16: 3802/3891, 5: 5309/5320, 6: 5422/5426, 7: 2898/2898, 12: 3319/3408, 13: 2380/2402, 14: 1923/1925, 15: 164/166. The pass rule is applied to the largest post-infusion table only.
- Freeze rule kept: only headers, barcodes and counts were printed; no CDR3, gene or clonotype value was printed or tabulated. Only chain counts (TRB 4,505, TRA 3,886) and the number of distinct `clonotype_ids` (5,027) of that one table were read, as counts.

## Formats

- **TCR tables** (inside `GSE125881_RAW.tar`, 16 members, 2,048,000 bytes, one `*-unpaired-clonotypes.csv.gz` per sample, 6,905 to 305,886 bytes; the spec's "about 7 MB" does not match). Columns (13, all read as text): `chain, clonotype_ids, consensus_ids, barcodes, contig_ids, cdr3s, cdr3_nts, v_genes, d_genes, j_genes, c_genes, frequency, proportion`. One row per clonotype x chain (TRA or TRB), not per cell, and not 10x contig annotations: there is no productive/full-length flag.
- **Expression matrix** `GSE125881_raw.expMatrix.csv.gz` (listed 109M, 114,573,707 bytes by HEAD): CSV, **genes x cells**. Its first line is an empty first field followed by 62,167 quoted cell IDs `"<16-nt>-<aggr ID>"`, aggr IDs 1 to 16. Rows are therefore genes, and Step 5's whole-matrix fetch was not needed. Cells per aggr ID: 1: 5515, 2: 6722, 3: 4920, 4: 1736, 5: 6117, 6: 6488, 7: 5023, 8: 1763, 9: 6271, 10: 4287, 11: 669, 12: 3484, 13: 2459, 14: 2382, 15: 197, 16: 4134.
- NHL-6 IP and NHL-7 IP were sequenced on the same chromium index; GEO says VDJ cell barcodes were used to tell them apart. Their tables are assigned by aggr ID (16 and 12), which the join above confirms.

## GSM table (GSM x patient x time point; the 16 VDJ GSMs are the GE GSM + 16)

| GE GSM | VDJ GSM | Patient | Time point as written | Stage label | Day | aggr ID |
|---|---|---|---|---|---|---|
| GSM3583861 | GSM3583877 | CLL-1 | IP | IP | 0 (infusion product) | 1 |
| GSM3583862 | GSM3583878 | CLL-1 | d21 | Early | 21 | 2 |
| GSM3583863 | GSM3583879 | CLL-1 | d38 | Contraction | 38 | 3 |
| GSM3583864 | GSM3583880 | CLL-1 | d112 | Late | 112 | 4 |
| GSM3583865 | GSM3583881 | CLL-2 | IP | IP | 0 (infusion product) | 8 |
| GSM3583866 | GSM3583882 | CLL-2 | d12 | Early | 12 | 9 |
| GSM3583867 | GSM3583883 | CLL-2 | d29 | Contraction | 29 | 10 |
| GSM3583868 | GSM3583884 | CLL-2 | d83 | Late | 83 | 11 |
| GSM3583869 | GSM3583885 | NHL-6 | IP | IP | 0 (infusion product) | 16 |
| GSM3583870 | GSM3583886 | NHL-6 | d12 | Early | 12 | 5 |
| GSM3583871 | GSM3583887 | NHL-6 | d29 | Contraction | 29 | 6 |
| GSM3583872 | GSM3583888 | NHL-6 | d102 | Late | 102 | 7 |
| GSM3583873 | GSM3583889 | NHL-7 | IP | IP | 0 (infusion product) | 12 |
| GSM3583874 | GSM3583890 | NHL-7 | d12 | Early | 12 | 13 |
| GSM3583875 | GSM3583891 | NHL-7 | d28 | Contraction | 28 | 14 |
| GSM3583876 | GSM3583892 | NHL-7 | d89 | Late | 89 | 15 |

"Day" is the number in the label; "IP" has no day number, and 0 is only a placeholder for "infusion product" (the product is not a post-infusion day). The "Stage label" is the GEO field "time point - 2". Four patients, four time points each.

## Licence lines

- Article: CC BY 4.0 (PMC6954177, Sheih et al. 2020, Nat Commun 11:219, doi:10.1038/s41467-019-13880-1).
- Data: public GEO accession GSE125881. GEO/NCBI terms were not read here: the GEO HTML pages sit behind a CAPTCHA, which is not bypassed. UNVERIFIED; cite the paper and the accession.
- Cao 2025 (CC BY-NC-ND 4.0, PMC12055983) was not used.

## Files and hashes (this machine, session scratchpad; not the dataset cache)

- Tar: `/private/tmp/claude-501/-Users-shreyan-lumenCell-V2/453b6f74-ae4a-4614-a97a-fbc583e58d2c/scratchpad/sheih-verify/GSE125881_RAW.tar`, 2,048,000 bytes, SHA-256 `769fa5b71b12b7efe433f51164d4370e0af00c20bebc647029cc4fbd9a943487`.
- Chosen member `GSM3583882_CLL-2_d12-unpaired-clonotypes.csv.gz` (305,886 bytes; read in memory with `tar -xO`, never extracted): SHA-256 `93f8ef95bbd990ee92b280bf0aa3922821c702bea3494898df3e0292dbcfa5db`; 1,832,439 bytes uncompressed.
- Matrix header (uncompressed first line, 1,323,122 bytes): `.../scratchpad/sheih-verify/header.txt`, SHA-256 `23955ade7cc347f34acc9668a537caab01f56a6126112abab5de529cfcf6315f`. Cell IDs with suffix `-9`: `.../scratchpad/cells_9.txt`.
- They are not in `~/Library/Caches/lumencell/datasets` (the fetch used a fresh scratchpad directory). The scratchpad may be cleaned; Tasks 10 and 11 fetch `GSE125881_RAW.tar` again (approved) and re-verify the SHA-256 above.

## Handling

Downloaded: `GSE125881_RAW.tar` into a fresh directory with `fetch_source` (certifi, size check, atomic rename). Its members were read in memory and never extracted to disk; the script read a member from a pipe. The matrix was streamed for its first line only (393,254 compressed bytes read, under the 16 MB cap).
