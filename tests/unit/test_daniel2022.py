from lumencell.bench.builders.daniel2022 import GEO_SUPPL, LICENCE, SOURCES

EXPECTED = (
    "GSE188666_scRNA_LCMV.mtx.gz",
    "GSE188666_scRNA_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_LCMV_features.tsv.gz",
    "GSE188666_scRNA_LCMV_metadata.tsv.gz",
    "GSE188666_scRNA_D5_LCMV.mtx.gz",
    "GSE188666_scRNA_D5_LCMV_barcodes.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_features.tsv.gz",
    "GSE188666_scRNA_D5_LCMV_metadata.txt.gz",
    "GSE188666_RAW.tar",
)


def test_sources_are_the_nine_geo_files() -> None:
    assert tuple(s.filename for s in SOURCES) == EXPECTED
    assert all(s.url == GEO_SUPPL + s.filename and s.sha256 and len(s.sha256) == 64 for s in SOURCES)
    assert 1.0e9 < sum(s.size for s in SOURCES) < 1.2e9 and not any("GSE188670" in s.url for s in SOURCES)


def test_the_licence_is_recorded() -> None:
    assert "GSE188666" in LICENCE and "10.1038/s41590-022-01337-5" in LICENCE
