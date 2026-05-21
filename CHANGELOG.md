# Changelog

## 2026-05-21

- Published Docker image `jinyuede163/scpipline-rstudio:singlecell-pipeline-20260521` (`sha256:f6814eed9ca019752e56034d8ec556d64ea31b5e1df582ff2dcec3edc993087c`).
- Updated the fastp/Cell Ranger workflow for Cell Ranger 10 by adding `create_bam` and `include_introns` configuration support.
- Replaced Cell Ranger `outs` symlink staging with explicit copies of filtered matrix MEX, filtered matrix HDF5, web summary HTML, and metrics summary CSV.
- Synced downstream single-cell analysis compatibility updates, including Seurat v5-related fixes and refreshed helper scripts.
- Added `docker/data/config_Example.ini` as a runnable configuration template.
