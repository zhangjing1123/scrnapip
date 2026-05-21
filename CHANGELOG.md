# Changelog

## 2026-05-21

- Published Docker image `jinyuede163/scpipline-rstudio:singlecell-pipeline-20260521` (`sha256:f6814eed9ca019752e56034d8ec556d64ea31b5e1df582ff2dcec3edc993087c`).
- Updated the fastp/Cell Ranger workflow for Cell Ranger 10 by adding `create_bam` and `include_introns` configuration support.
- Replaced Cell Ranger `outs` symlink staging with explicit copies of filtered matrix MEX, filtered matrix HDF5, web summary HTML, and metrics summary CSV.
- Synced downstream single-cell analysis compatibility updates, including Seurat v5-related fixes and refreshed helper scripts.
- Documented the workflow change from Cerebro export to Loupe Browser `.cloupe` export via `loupeR`.
- Documented the replacement of the old CytoTRACE step with CytoTRACE2 outputs and configuration.
- Added `docker/data/config_Example.ini` as a runnable configuration template.
- Confirmed the bundled `fastp` binary is unchanged from the maintained pipeline and remains fastp 0.23.2.
- Synced additional downstream fixes: DoubletFinder runs per sample before integration, Seurat v5 layer joins/fallbacks are used where needed, Monocle 2 compatibility patches are included, CellChat optional `projectData()` handling is guarded, clusterProfiler database lookup is stabilized, and step12 Circos output is restored.
