"""Daniel et al. 2022 (chronic LCMV, paired scRNA and TCR): the GEO files."""

from lumencell.bench.download import Source

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE188nnn/GSE188666/suppl/"

_FILES: tuple[tuple[str, int], ...] = (
    ("GSE188666_scRNA_LCMV.mtx.gz", 602240886),
    ("GSE188666_scRNA_LCMV_barcodes.tsv.gz", 487511),
    ("GSE188666_scRNA_LCMV_features.tsv.gz", 113716),
    ("GSE188666_scRNA_LCMV_metadata.tsv.gz", 4846117),
    ("GSE188666_scRNA_D5_LCMV.mtx.gz", 544512106),
    ("GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz", 389725),
    ("GSE188666_scRNA_D5_LCMV_features.tsv.gz", 113719),
    ("GSE188666_scRNA_D5_LCMV_metadata.txt.gz", 3976768),
    ("GSE188666_RAW.tar", 12113920),
)

# SHA-256 as recorded on the first download (2026-11-19). The two mtx sizes are 574 and 519 MiB.
_SHA256 = {
    "GSE188666_scRNA_LCMV.mtx.gz": ("be4a6e2825d9cdffbb4bebb094e43be00ed388faa11a5bac95c2418bdd576ed6"),
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz": (
        "c75bc8d63c793b4c74215d46d40d243118b5fb8f4e1aff248af0451ee15196fe"
    ),
    "GSE188666_scRNA_LCMV_features.tsv.gz": (
        "dee1b30b23a9ceae7117d71b894ec7e4d5441fc934fbcb242493ef57b2e06e2d"
    ),
    "GSE188666_scRNA_LCMV_metadata.tsv.gz": (
        "0594a0463b0d667307736fa45453d9aae9ed761ce105d53b3a24c30dffa46393"
    ),
    "GSE188666_scRNA_D5_LCMV.mtx.gz": ("14bc8bf47f68d11b27dbabfc8ec6613385f58a83578a0a81b7187817602657b3"),
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz": (
        "3ec9ca5d7650e5ab428121161c7533d1da42d72ca1187dd9941e8de0207151c8"
    ),
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz": (
        "e3367683bbf603ef45124593bcc9a2868ffe393e422df1a5a035c24e7e76a576"
    ),
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz": (
        "9f8943dea624e8f56c894421eaad3b723b9e2698a3f8d148122fa516bab7c4ec"
    ),
    "GSE188666_RAW.tar": ("db4932008b1b9e603c2f076aa7cf60302f38d8b66108188f3827a121a676c9a9"),
}

SOURCES: tuple[Source, ...] = tuple(Source(GEO_SUPPL + n, n, size, _SHA256[n]) for n, size in _FILES)

LICENCE = (
    "Data: GEO GSE188666 (SuperSeries GSE188670); "
    "UNVERIFIED: https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html is behind a CAPTCHA, "
    "not bypassed (read attempt 2026-11-19). "
    "Paper: Daniel et al. 2022, Nat Immunol 23:1614-1627, doi:10.1038/s41590-022-01337-5; "
    "UNVERIFIED: https://pmc.ncbi.nlm.nih.gov/articles/PMC11225711/ is behind a CAPTCHA, "
    "not bypassed (read attempt 2026-11-19)."
)
