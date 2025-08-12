from __future__ import annotations
from pathlib import Path
import re
import anndata as ad

inputs = list(snakemake.input)
out_path = Path(snakemake.output[0])
batch_key = snakemake.params.get("batch_key", "sample")

adatas = []
keys = []
for p in inputs:
    p = Path(p)
    # infer sample ID from "results/qc/{sample}.filtered.h5ad"
    m = re.match(r"(.+)\.filtered\.h5ad$", p.name)
    sample = m.group(1) if m else p.stem
    keys.append(sample)
    A = ad.read_h5ad(p)
    # prefer concat to assign batch column; but ensure column exists if user wants it explicit
    A.obs[batch_key] = sample
    adatas.append(A)

combined = ad.concat(adatas, label=batch_key, keys=keys, join="outer", fill_value=0)
out_path.parent.mkdir(parents=True, exist_ok=True)
combined.write_h5ad(out_path, compression="gzip")

# from pathlib import Path
# import anndata as ad

# out = Path(snakemake.output[0])  # results/integration/concat.h5ad
# out.parent.mkdir(parents=True, exist_ok=True)

# adatas = []
# for p in snakemake.input:  # per-sample filtered H5ADs
#     a = ad.read_h5ad(p)
#     # sample name from filename (e.g. S01.filtered.h5ad -> S01)
#     sample = Path(p).stem.split(".")[0].split("_")[0]
#     a.obs["sample"] = sample
#     a.var_names_make_unique()
#     adatas.append(a)

# cat = ad.concat(adatas, join="inner", label="sample", merge="unique", index_unique=None)

# # also add a configurable batch key if you want something other than 'sample'
# bk = snakemake.params.get("batch_key", "sample")
# if bk != "sample":
#     cat.obs[bk] = cat.obs["sample"].values

# cat.write_h5ad(out, compression="gzip")