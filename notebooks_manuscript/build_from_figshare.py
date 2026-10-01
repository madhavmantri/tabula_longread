#!/usr/bin/env python3
"""Rebuild the derived objects the analysis notebooks need from the Figshare deposit.

Run from the repository root, after placing the Figshare files in pacbio/h5ads/ and csvs/:

    python notebooks_manuscript/build_from_figshare.py prepare
    # ... run 01_atlas_construction/01_09_compute_isoform_fraction.ipynb ...
    python notebooks_manuscript/build_from_figshare.py cdkn2a

`prepare`
  - Writes the annotated-only level, which is the ensemblid level restricted to features
    carrying an Ensembl transcript identifier (ENST...). This is the same feature set and
    cell set as the object used in the manuscript.
  - Links the deposited classification tables to the file names the notebooks read.
    The deposited tables hold one row for every feature in the pbid object (the rows the
    notebooks join on), with SQANTI3's ORF_length renamed ORF_genomic_span because it
    reports the genomic span of the CDS, not its length.
  - Writes the per-cell popV predictions from the pbid object's obs.
  - Derives the ensemblid-level SQANTI3 table from the pbid-level one. This is close to,
    but not identical to, the manuscript table (see the comment in prepare()).

`cdkn2a`
  Writes the *_xgen_cdkn2a.h5ad objects read by 06_senescence/: each input object restricted
  to the cells of the four capture donors, with the CDKN2A capture calls from
  xgen_cdkn2a_per_cell_pigeon_raw.csv attached as obs columns. The raw-count columns are
  identical to those used in the manuscript. The distinct-UMI columns (xgen_*_umi) are not
  deposited and are not written; only the recovery comparison in 06_03 uses them.
"""
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

H5AD = Path("pacbio/h5ads")
CSV = Path("csvs")
STEM = "all_samples_pacbio_recollapsed_raw_counts_{level}_bc_anndata_preprocessed"

P16 = ["ENST00000304494", "ENST00000494262", "ENST00000498628",
       "ENST00000498124", "ENST00000578845", "ENST00000380151"]
P14ARF = ["ENST00000579755", "ENST00000530628"]

# Objects the senescence notebooks load, as the input stem each is built from.
CDKN2A_INPUTS = [
    STEM.format(level="genes"),
    STEM.format(level="pbids"),
    STEM.format(level="ensemblids"),
    STEM.format(level="ensemblids_annotatedonly"),
    STEM.format(level="pbids") + "_with_isofrac",
    STEM.format(level="ensemblids") + "_with_isofrac",
    STEM.format(level="ensemblids_annotatedonly") + "_with_isofrac",
]


def _exists(path):
    """Never overwrite: an existing file may be the original, not a rebuild."""
    if path.exists():
        print(f"skip {path.name} (already exists)")
        return True
    return False


def prepare():
    src = H5AD / f"{STEM.format(level='ensemblids')}.h5ad"
    dst = H5AD / f"{STEM.format(level='ensemblids_annotatedonly')}.h5ad"
    if not _exists(dst):
        adata = ad.read_h5ad(src)
        keep = adata.var_names.str.startswith("ENST")
        adata[:, keep].copy().write_h5ad(dst)
        print(f"wrote {dst.name}: {adata.n_obs:,} cells x {int(keep.sum()):,} features")
        del adata

    for tool in ("pigeon", "sqanti3"):
        deposited = CSV / f"all_samples_pbids_feature_metadata_{tool}.csv"
        expected = CSV / f"all_samples_recollapsed_{tool}_classification.csv"
        if not expected.exists() and not expected.is_symlink():
            expected.symlink_to(deposited.name)
        print(f"{expected.name} -> {deposited.name}")

    pbids = ad.read_h5ad(H5AD / f"{STEM.format(level='pbids')}.h5ad", backed="r")
    var, obs = pbids.var[["transcript_id"]].copy(), pbids.obs.copy()
    pbids.file.close()

    # popV labels are stored on the objects; write them out in the table layout 03_02 reads.
    dst = CSV / "all_samples_recollapsed_popv_predictions_df.csv"
    if not _exists(dst):
        obs[["popv_prediction", "popv_prediction_score"]].to_csv(dst)
        print(f"wrote {dst.name}: {len(obs):,} cells")

    # Ensemblid-level SQANTI3 table: one row per transcript_id (ENST where the structure
    # matched an annotated transcript, else its PB id), keeping the first pbid row of each,
    # with `isoform` set to the transcript_id so add_classification() joins on the
    # ensemblid object. The manuscript chose that first row among all pre-filter structures;
    # only the atlas pbids are deposited, so for about 16% of ENST features a different
    # structure represents the transcript and its annotation can differ. PB features are
    # identical to the manuscript table.
    dst = CSV / "all_samples_recollapsed_sqanti3_classification_ensemblids.csv"
    if not _exists(dst):
        cls = pd.read_csv(CSV / "all_samples_pbids_feature_metadata_sqanti3.csv",
                          dtype=str, keep_default_na=False)
        cls["transcript_id"] = cls["isoform"].map(var["transcript_id"])
        cls = (cls.dropna(subset=["transcript_id"])
                  .drop_duplicates(subset="transcript_id")
                  .rename(columns={"isoform": "molecule_id"}))
        cls["isoform"] = cls.pop("transcript_id")
        cls.to_csv(dst, index=False)
        print(f"wrote {dst.name}: {len(cls):,} transcript_ids")


def cdkn2a():
    calls = pd.read_csv(CSV / "xgen_cdkn2a_per_cell_pigeon_raw.csv", index_col="cell_id")
    p16, p14 = calls[P16].sum(axis=1), calls[P14ARF].sum(axis=1)
    obs_cols = pd.DataFrame({
        "xgen_matched": calls["xgen_matched"],
        "xgen_cdkn2a_raw": calls["cdkn2a_total"],
        "xgen_p16_raw": p16,
        "xgen_p14ARF_raw": p14,
        "xgen_other_cdkn2a_raw": calls["cdkn2a_total"] - p16 - p14,
        "xgen_cdkn2a_raw_pos": calls["cdkn2a_pos"],
        "xgen_cdkn2a_raw_status": calls["cdkn2a_status"],
    }, index=calls.index)

    for stem in CDKN2A_INPUTS:
        src = H5AD / f"{stem}.h5ad"
        if not src.exists():
            print(f"skip {src.name} (not found; run 01_09 first for *_with_isofrac)")
            continue
        dst = H5AD / f"{stem}_xgen_cdkn2a.h5ad"
        if _exists(dst):
            continue
        adata = ad.read_h5ad(src)
        missing = obs_cols.index.difference(adata.obs_names)
        if len(missing):
            raise ValueError(f"{src.name}: {len(missing):,} capture cells absent from the object")
        # Keep the input object's cell order, as the manuscript objects do.
        adata = adata[adata.obs_names.isin(obs_cols.index)].copy()
        for col in obs_cols:
            adata.obs[col] = obs_cols.loc[adata.obs_names, col].to_numpy()
        adata.write_h5ad(dst)
        print(f"wrote {dst.name}: {adata.n_obs:,} cells")


if __name__ == "__main__":
    steps = {"prepare": prepare, "cdkn2a": cdkn2a}
    if len(sys.argv) != 2 or sys.argv[1] not in steps:
        sys.exit(__doc__)
    steps[sys.argv[1]]()
