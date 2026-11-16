# fap: fast AlphaPulldown

fap builds, submits and post-processes AlphaPulldown and AlphaFold 3 jobs that
pair GPCRs with candidate ligands on the Sherlock cluster. It was split out of
[ligandFinder](https://github.com/kbrulois/ligandFinder) on 2026-11-13, with
the history of the moved files.

- **fap depends on ligandFinder.** `library(fap)` attaches ligandFinder too, so
  a script gets everything `library(ligandFinder)` used to give it. fap's R code
  imports seven ligandFinder names: `get_db_path`, `map_table`,
  `factor_to_uniprotFeature`, `html_slide_show`, `distance`, `AA_colors` and
  `extra_color`. It also reads ligandFinder's `id_mapping.rds`, `gpcr_list.rds`,
  `ligand_list.rds` and GPCRdb pairings with
  `system.file(..., package = "ligandFinder")`.
- **ligandFinder does not depend on fap.** What crosses back is data: the pocket
  labels (`IC_dist` in ligandFinder's `knowns.rds`) came from these models.

## On the cluster

The cluster setup lives in one place, `inst/cluster/sherlock_env.sh`. It
exports `R_LIBS`, the group R library that holds fap and ligandFinder. Install
both packages into that library, then source the file from `~/.bashrc` once:

```bash
echo 'source /home/groups/ebutcher/programs/pipeline/R_libs4.1/fap/cluster/sherlock_env.sh' >> ~/.bashrc
```

Install or update both packages there, from Sherlock (ligandFinder first; fap depends on it):

```r
lib <- "/home/groups/ebutcher/programs/pipeline/R_libs4.1"
remotes::install_github("kbrulois/ligandFinder", lib = lib)
remotes::install_github("kbrulois/fap", lib = lib)
```

After that, scripts need only `library(fap)`. Shell lines reach installed
scripts as `$R_LIBS/fap/scripts/<folder>/<script>` or with
`system.file("scripts/<folder>/<script>", package = "fap")`. sbatch jobs inherit
`R_LIBS` from the shell that submits them.

## Where things are

| path | what |
|---|---|
| `R/afpd_*.R` | AlphaPulldown: job tables and submission (`afpd_submit_model_jobs.R`), features (`afpd_features.R`), AF3-style name parsing and renaming (`afpd_af3_name_parsing.R`), metric extraction (`afpd_metric_processing.R`) and assembly (`afpd_metric_assembly.R`) |
| `R/af3_metric_processing.R` | AlphaFold 3 output: file-name cleanup and raw metric import |
| `R/metric_retrieval.R` | reading extracted metrics and contacts back (`get_metrics`, `get_contacts`, `get_known_pairs`) |
| `R/initial_guess_pose.R` | initial-guess poses for AlphaPulldown |
| `R/af_db_access.R`, `R/cleave_peptide.R` | the local AlphaFold DB cache; cutting a peptide out of a predicted structure |
| `R/coev_*.R` | coevolution between receptor and ligand MSAs, and its plots |
| `R/lists_utils.R` | GPCR lists: Ballesteros-Weinstein numbering (`summarize_bw`, `getBallesterosFromGPCRDB`, ...) |
| `R/utils.R` | AlphaFold and GPCR constants: residue-pair classes, BW segments, file-type codes |
| `R/sysdata.rda` | `gpcr_list`, internal |
| `data/Chemokine_Pairs.rds` | known chemokine pairs |
| `inst/extdata/` | receptor and ligand lists that only fap uses (GPCR selections, FASTAs, peptide lists) |
| `inst/cluster/sherlock_env.sh` | the cluster setup (above) |
| `exec/` | command-line entry points: `afpd_submit_jobs_cli.R`, `afpd_rename_files_cli.R`, `af3_rename_files_cli.R`, `generate_spoc_json.py` |
| `demo/` | renaming demos |

### `inst/scripts/`, by pipeline step

| folder | what the scripts do |
|---|---|
| `jobs/` | build AlphaPulldown job tables (`afpd_create_job_input_dir.R`, `afpd_create_job_input_big_rec.R`), fix and restart jobs, submit them (`run_jobs.sh` -> `exec/afpd_submit_jobs_cli.R`) |
| `jobs/runs/` | one script per submitted batch: what was asked for, and when (CXCL14, GDF5/BRINP3, top-50 dibasic, uniprot peptides, ...) |
| `features/` | MSA and template features (`afpd_create_features.sh`, an sbatch array) and the initial-guess route (`alphapulldown_initial_guess.sh`, `ap_initial_guess_*.R`, the AlphaPulldown patch) |
| `metrics/` | raw output -> metrics: rename and gate (`afpd_prepare_metric_jobs.R`, `afpd_write_side.R`), the sbatch array (`extract_metrics.sh` -> `afpd_extract_metrics_worker.R`), assembly (`afpd_assemble_metrics.R`, `afpd_assemble_known_pairs.R`), ColabFold import, `parse_pkl.py`, `voronota.sh`, a test on one model |
| `contacts/` | contact representations of the complexes and the models trained on them (R reports, Python xgboost / feed-forward / attention models, SHAP, UMAP) |
| `af3/` | AlphaFold 3 runs: renaming, dataset checks, contact report |
| `files/` | housekeeping on the cluster: unnest / tar / download model directories, renaming, name-log checks |
| `structures/` | AlphaFold DB and PDB downloads (`af_db_bulk_download.sh` builds the cache `R/af_db_access.R` reads) |
| `receptor_lists/` | building and refining the GPCR and ligand lists |
| `coev/` | coevolution: fetch MSAs, run the analysis, reports |

### The main pipeline

1. **Jobs:** `jobs/afpd_create_job_input_dir.R` writes the job table, and
   `jobs/run_jobs.sh` submits it through `submit_model_jobs()`.
2. **Features:** `features/afpd_create_features.sh`, or the initial-guess
   route through `features/alphapulldown_initial_guess.sh`.
3. **Metrics:** `metrics/afpd_write_side.R` renames the raw directories,
   checks them and queues `metrics/extract_metrics.sh`. Each array task runs
   `afpd_extract_metrics_worker.R` -> `do_metrics()`, and
   `metrics/afpd_assemble_metrics.R` collects the results.
4. **Analysis:** `contacts/`, with `get_metrics()` / `get_contacts()` to read
   the results back.

## Known issues

- `af3/af3_contact_rep.R` and `metrics/import_colabfold_metrics.R` call
  `do_metrics(residue_data = ...)`. The function has taken `res_dat` since
  2026-04-02, so both scripts stop with an "unused argument" error.
- `License:` in DESCRIPTION is still the RStudio placeholder, as in ligandFinder.
- Three scripts do not parse, and did not in ligandFinder either. They are
  unfinished scratch: `files/afpd_unnest_dirs_tar.R`,
  `receptor_lists/refine_gpcr_list.R` and `structures/fetch_pdb_files.R`.
