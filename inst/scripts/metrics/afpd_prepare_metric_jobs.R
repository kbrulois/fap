#!/usr/bin/env Rscript
## ---------------------------------------------------------------------------
## afpd_prepare_metric_jobs.R
##
## Raw AlphaPulldown output -> renamed directories -> submitted metric jobs,
## for the three CXCL14 / VTM2A trees.  The read-only pre-flight
## (rename_dryrun.R) is assumed to have passed; this is the write side.
##
## Three steps per tree.  Run them one at a time, in order, and read the
## printout between them:
##
##   Rscript afpd_prepare_metric_jobs.R A survey           # read-only
##   Rscript afpd_prepare_metric_jobs.R A rename           # renames dirs+files
##   Rscript afpd_prepare_metric_jobs.R A queue            # writes group files
##   Rscript afpd_prepare_metric_jobs.R A queue --submit   # ... and sbatches
##
## Every step re-derives its state from the filesystem and keeps no memory of
## earlier passes, so each is safely re-runnable.  When the --mem=64G restart
## lands new directories in trees A and C, run `rename` and then `queue`
## again: only the new directories move, nothing already renamed is touched,
## and nothing already extracted is re-extracted.
##
## ALLOCATION.  `rename` and `queue` refuse to start without SLURM_JOB_ID,
## because both are too heavy for a login node:
##
##   salloc --partition=normal --cpus-per-task=8 --mem=64G --time=4:00:00
##
## `rename` reads and rewrites the 43M-row random_codes.csv (a few GB
## resident).  `queue` materialises every metrics_v2.csv and every
## metrics_v2c.rds into list columns, which is why it runs in blocks of
## QUEUE_CHUNK directories and drops the list columns in the same expression
## that computes the flags.  `survey` only lists directories and reads
## ranking_debug.json, so it is fine on a login node.
##
## Do not run two trees at once: do_renaming() rewrites random_codes.csv in
## place, and concurrent passes would hand out the same codes twice.
##
## QUARANTINE.  A ranking_debug.json that is unreadable, or that lacks "order"
## / "iptm" / "iptm+ptm", kills do_renaming() in parse_afpd_files() -- after
## stage 3 has already renamed every directory in the batch.  Such directories
## are therefore excluded from the batch BEFORE do_renaming() is called and
## listed in a quarantine report.  They are excluded by name rather than moved
## aside: moving them would stop AlphaPulldown resuming that pair in place, and
## nothing except this script reads the raw tree.
## ---------------------------------------------------------------------------

.libPaths("/home/groups/ebutcher/programs/pipeline/R_libs4.1")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tidyr)
  library(ligandFinder)
})

## --- per-tree parameters (requirement: one script, no fork for run B) ------
##
## run_name       goes into the new file names (do_renaming) and must equal
##                basename(run_id_dir), see below.
## run_id_dir     submit_metric_jobs() takes run_id = fs::path_file(run_dir)
##                and writes it as the 4th column of the group file; the worker
##                passes it to do_metrics(run_name=).  Tree B's model
##                directories sit in <run>/models, whose basename is "models",
##                so B points one level up and records "cxl14_75Y" like the
##                other two.  Nothing else reads this column.
## expected_models  the ABSOLUTE number of predictions per directory: 5 for the
##                standard driver (--num_predictions_per_model=1), 25 for run B
##                (--num_predictions_per_model=5).  This is the number both the
##                rename gate and the queue gate test against, so a short
##                prediction set can never be mistaken for a finished one.

TREES <- list(
  A = list(run_dir         = "/scratch/groups/ebutcher/deorphan/models/cxl14_x75",
           run_name        = "cxl14_x75",
           run_id_dir      = "/scratch/groups/ebutcher/deorphan/models/cxl14_x75",
           expected_models = 5L),
  B = list(run_dir         = "/oak/stanford/groups/ebutcher/kevin/cxl14_75Y/models",
           run_name        = "cxl14_75Y",
           run_id_dir      = "/oak/stanford/groups/ebutcher/kevin/cxl14_75Y",
           expected_models = 25L),
  C = list(run_dir         = "/scratch/groups/ebutcher/deorphan/models/vtm2a_157",
           run_name        = "vtm2a_157",
           run_id_dir      = "/scratch/groups/ebutcher/deorphan/models/vtm2a_157",
           expected_models = 5L)
)

RC_FILE     <- "/home/groups/ebutcher/kevin/ligandFinder/random_codes.csv"
REPORT_ROOT <- "/scratch/groups/ebutcher/deorphan/ligandFinder/prep_reports"
JOB_ROOT    <- "/scratch/groups/ebutcher/deorphan/ligandFinder/metric_jobs"
GROUP_SIZE  <- 100L
QUEUE_CHUNK <- 200L

## Shape of a renamed directory: exactly one "_", "h" + receptor, "h" + ligand
## + x<first>x<last>.  Used both to recognise renamed directories and to assert
## the rename targets before renaming anything.  Entry Names are upper case and
## digits, so this also rejects a target whose receptor failed to resolve (the
## name would come out as "h_hCXL14x100x109") and any name carrying a lower
## case "x", which would make parse_ranges() cut in the wrong place.
TARGET_RX <- "^h[A-Z0-9]+_h[A-Z0-9]+x\\d+x\\d+$"

## The naming arguments.  The same list is used for the dry computation of the
## targets and for the do_renaming() call, so what is asserted is what runs.
RENAME_ARGS <- list(delim_proteins  = "_",
                    delim_ranges    = "x",
                    delim_start_end = "x",
                    p1_prefix       = "h",
                    p1_suffix       = NA,
                    p2_prefix       = "h",
                    exclude_p1_range = TRUE)

CODE_ARGS <- list(site      = "SU",
                  submitter = "KB",
                  algorithm = "AF2v3",
                  random_seed = 42,
                  rc_file   = RC_FILE)

STAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")

## --- helpers --------------------------------------------------------------

write_report <- function(key, name, dat) {
  d <- fs::path(REPORT_ROOT, key)
  fs::dir_create(d)
  f <- fs::path(d, paste0(STAMP, "_", name, ".tsv"))
  data.table::fwrite(as.data.frame(dat), f, sep = "\t")
  message("report: ", f, "  (", nrow(dat), " rows)")
  invisible(f)
}

require_allocation <- function(step) {
  if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) {
    stop("'", step, "' must not run on a login node (it would be OOM-killed).\n",
         "  salloc --partition=normal --cpus-per-task=8 --mem=64G --time=4:00:00\n",
         "  ml R/4.1.2", call. = FALSE)
  }
}

## future::plan() does not evaluate `workers` where it is written -- a local
## variable is invisible to it -- so the count is substituted into the call.
set_future_plan <- function(workers) {
  eval(bquote(future::plan(strategy = future::multicore(workers = .(as.integer(workers))))))
}

count_files <- function(dir_path, pattern) {
  vapply(dir_path, function(d) length(list.files(d, pattern = pattern)),
         integer(1), USE.NAMES = FALSE)
}

## ranking_debug.json, read exactly the way the pipeline reads it:
## jsonlite::fromJSON() as in parse_afpd_files(), and the three fields
## parse_afpd_files() ("order") and import_raw_metrics() ("iptm", "iptm+ptm")
## require.  Returns one row per directory.
rd_status <- function(dir_path) {
  f <- fs::path(dir_path, "ranking_debug.json")
  na <- list(n_iptm = NA_integer_, n_ipptm = NA_integer_, n_order = NA_integer_)
  if (!file.exists(f)) return(tibble(rd = "absent", !!!na))
  j <- tryCatch(jsonlite::fromJSON(f), error = function(e) NULL)
  if (is.null(j) || !is.list(j)) return(tibble(rd = "unparseable", !!!na))
  miss <- setdiff(c("iptm", "iptm+ptm", "order"), names(j))
  if (length(miss)) {
    return(tibble(rd = paste0("missing_", paste(miss, collapse = "+")), !!!na))
  }
  tibble(rd      = "ok",
         n_iptm  = length(j[["iptm"]]),
         n_ipptm = length(j[["iptm+ptm"]]),
         n_order = length(j[["order"]]))
}

## --- step 1: survey (read-only, no filesystem effect) ---------------------

survey_tree <- function(tree) {
  E <- tree[["expected_models"]]
  dirs <- fs::dir_ls(tree[["run_dir"]], type = "directory")
  if (length(dirs) == 0) stop("no directories in ", tree[["run_dir"]], call. = FALSE)

  out <- tibble(afpd_dir_name = as.character(fs::path_file(dirs)),
                afpd_dir      = as.character(dirs)) %>%
    mutate(shape = case_when(
      stringr::str_detect(afpd_dir_name, "_and_")    ~ "raw",
      stringr::str_detect(afpd_dir_name, TARGET_RX)  ~ "renamed",
      TRUE                                           ~ "unknown"))

  out <- bind_cols(out, map_dfr(out[["afpd_dir"]], rd_status))

  out %>%
    mutate(stage = case_when(
      shape == "unknown"                      ~ "unknown_name",
      rd == "absent"                          ~ "waiting_af",
      rd != "ok"                              ~ "quarantine_json",
      n_iptm != E | n_ipptm != E | n_order != E ~ "short_models",
      shape == "raw"                          ~ "rename_ready",
      TRUE                                    ~ "renamed"))
}

step_survey <- function(tree, key) {
  s <- survey_tree(tree)
  message("\n", key, "  ", tree[["run_dir"]],
          "   expected models per directory: ", tree[["expected_models"]])
  print(count(s, shape, stage, rd))
  write_report(key, "survey", s %>% select(-afpd_dir))
  invisible(s)
}

## --- step 2: rename -------------------------------------------------------

step_rename <- function(tree, key) {
  require_allocation("rename")
  E <- tree[["expected_models"]]

  s <- survey_tree(tree)
  print(count(s, shape, stage))
  write_report(key, "survey_before_rename", s %>% select(-afpd_dir))

  ## Quarantine (excluded by name, never handed to do_renaming) and the two
  ## other kinds of raw directory that are not rename-ready.  waiting_af and
  ## short_models are left raw on purpose: AlphaPulldown can still resume
  ## those pairs in place, and the next pass picks them up when they are whole.
  for (st in c("quarantine_json", "short_models", "unknown_name")) {
    bad <- s %>% filter(stage == st)
    if (nrow(bad)) {
      message("!! ", nrow(bad), " directories excluded as ", st)
      write_report(key, st, bad %>% select(-afpd_dir))
    }
  }

  elig <- s %>% filter(shape == "raw", stage == "rename_ready")
  message(nrow(elig), " raw directories are rename-ready; ",
          sum(s[["stage"]] == "waiting_af"), " still waiting on AlphaFold")
  if (nrow(elig) == 0) {
    message("nothing to rename")
    return(invisible(NULL))
  }

  ## The target names, computed with the pure exported functions -- zero
  ## filesystem effect -- from the same arguments do_renaming() will use.
  plan <- do.call(make_new_dirname,
                  c(list(input = parse_dirname(pairing_dir = elig[["afpd_dir_name"]])),
                    RENAME_ARGS))

  ## ---- assertions in front of the irreversible rename ------------------
  stopifnot(nrow(plan) == nrow(elig),
            setequal(plan[["afpd_dir_name"]], elig[["afpd_dir_name"]]))

  bad_shape <- plan %>%
    filter(is.na(new_dir_name) | !stringr::str_detect(new_dir_name, TARGET_RX))
  dup <- plan %>% count(new_dir_name) %>% filter(n > 1)
  ## file.rename() onto a non-empty target returns FALSE and warns; onto an
  ## EMPTY target it succeeds and replaces it.  Neither may happen.
  clash <- plan %>%
    filter(file.exists(fs::path(tree[["run_dir"]], new_dir_name)))

  if (nrow(bad_shape) || nrow(dup) || nrow(clash)) {
    if (nrow(bad_shape)) write_report(key, "FAIL_bad_target_name",
                                      bad_shape %>% select(afpd_dir_name, new_dir_name))
    if (nrow(dup))       write_report(key, "FAIL_duplicate_target", dup)
    if (nrow(clash))     write_report(key, "FAIL_target_exists",
                                      clash %>% select(afpd_dir_name, new_dir_name))
    stop("refusing to rename: ", nrow(bad_shape), " malformed target names, ",
         nrow(dup), " targets claimed by more than one directory, ",
         nrow(clash), " targets already on disk -- see the FAIL_* reports",
         call. = FALSE)
  }

  stopifnot(all(fs::dir_exists(fs::path(tree[["run_dir"]], plan[["afpd_dir_name"]]))),
            file.access(tree[["run_dir"]], 2) == 0)
  ## get_codes() freads random_codes.csv and fwrites it back from inside
  ## do_renaming(); if that write failed the run would abort with the
  ## directories renamed and the files not.
  stopifnot(file.exists(RC_FILE), file.access(RC_FILE, 2) == 0)

  ## Provenance audit, on disk BEFORE anything moves.
  write_report(key, "rename_plan",
               plan %>% select(afpd_dir_name, new_dir_name, model_name, model_id))

  do.call(do_renaming,
          c(list(run_dir     = tree[["run_dir"]],
                 run_name    = tree[["run_name"]],
                 pairing_dir = elig[["afpd_dir_name"]],
                 afpd_raw    = FALSE),
            RENAME_ARGS, CODE_ARGS))

  ## ---- verify what actually happened -----------------------------------
  ver <- plan %>%
    transmute(afpd_dir_name, new_dir_name,
              target = as.character(fs::path(tree[["run_dir"]], new_dir_name))) %>%
    mutate(target_exists = fs::dir_exists(target),
           source_gone   = !fs::dir_exists(fs::path(tree[["run_dir"]], afpd_dir_name)),
           n_ark         = count_files(target, "_ark_.*\\.pdb$"),
           has_file_log  = file.exists(fs::path(target, "file_name_log.csv"))) %>%
    select(-target)

  print(ver %>% count(target_exists, source_gone, has_file_log, n_ark))
  write_report(key, "rename_verify", ver)
  if (any(!ver[["target_exists"]] | !ver[["source_gone"]] | ver[["n_ark"]] != E)) {
    message("!! some directories did not land as planned -- read the ",
            "rename_verify report before running 'queue'")
  }
  invisible(ver)
}

## --- step 3: queue --------------------------------------------------------

## One block of directories: the package's own checks, then the flags, then
## the list columns are gone (metrics and contacts are the multi-GB ones).
check_block <- function(block, bw_rows, E) {
  afpd_check_metrics(block) %>%
    afpd_check_contacts() %>%
    mutate(metrics_rows = map_int(metrics, function(m) if (is.null(m)) 0L else nrow(m)),
           has_depth    = map_lgl(metrics, function(m) !is.null(m) && "depth" %in% colnames(m)),
           bw_rows      = unname(bw_rows[p1_name])) %>%
    select(afpd_dir_name, afpd_dir, run_dir, p1_name, p2_name,
           n_iptm, num_ark, num_con, num_pae,
           has_json_debug, has_metrics_csv, has_contact_rds, contacts_good,
           metrics_rows, has_depth, bw_rows)
}

step_queue <- function(tree, key, submit = FALSE) {
  require_allocation("queue")
  E <- tree[["expected_models"]]

  s <- survey_tree(tree)
  print(count(s, shape, stage))
  write_report(key, "survey_before_queue", s %>% select(-afpd_dir))

  ren <- s %>% filter(shape == "renamed")
  if (nrow(ren) == 0) {
    message("no renamed directories in ", tree[["run_dir"]], " -- run 'rename' first")
    return(invisible(NULL))
  }

  ncores <- max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4")))
  set_future_plan(max(1L, ncores %/% 2))
  options(future.globals.maxSize = 32 * 1024^3)

  ## zero-row "bw: full_table" (455 of 845 rows) and receptors absent from
  ## gpcr_list can never produce a depth column: process_metrics() derives the
  ## receptor sequence from that table, so seq_match2 is "different", nothing
  ## reaches compute_RLdists(), and do_metrics() writes no depth.  Structurally
  ## incapable, not unattempted.
  gl <- getFromNamespace("gpcr_list", "ligandFinder")
  bw_rows <- setNames(vapply(gl[["bw: full_table"]],
                             function(z) as.integer(NROW(z)), integer(1)),
                      gl[["uniprot_name"]])

  runs <- ren %>%
    transmute(afpd_dir_name, afpd_dir, n_iptm,
              run_dir = tree[["run_id_dir"]])

  blocks <- split(runs, ceiling(seq_len(nrow(runs)) / QUEUE_CHUNK))
  flags <- map_dfr(blocks, check_block, bw_rows = bw_rows, E = E)
  rm(blocks, runs); gc()

  flags <- flags %>%
    mutate(
      ## ABSOLUTE expected count, cross-checked against ranking_debug.json's
      ## iptm, not merely num_ark == num_con == num_pae.
      af_ready = has_json_debug & !is.na(n_iptm) &
                 n_iptm == E & num_ark == E & num_con == E & num_pae == E,
      done     = has_metrics_csv & metrics_rows == E & has_depth & contacts_good,
      ## do_metrics() writes metrics_v2.csv only on success, so has_metrics_csv
      ## means an extraction ran to the end: if it still has no depth column,
      ## or fewer rows than models, repeating it gives the same answer.
      state = case_when(
        done                                        ~ "done",
        num_ark == 0                                ~ "files_not_renamed",
        is.na(bw_rows) | bw_rows == 0               ~ "incapable_no_bw",
        has_metrics_csv & metrics_rows == E & !has_depth ~ "incapable_no_depth",
        has_metrics_csv & metrics_rows != E          ~ "incapable_short_metrics",
        ## a renamed directory is one AlphaFold already finished, so a model
        ## count below expected_models is a short or damaged prediction set,
        ## not work in progress: reported, never queued, never called done.
        !af_ready                                   ~ "incomplete_models",
        has_metrics_csv                             ~ "queue_retry",
        TRUE                                        ~ "queue_fresh"))

  message("\n", key, "  ", nrow(flags), " renamed directories")
  print(count(flags, state))
  write_report(key, "queue_state", flags)

  q <- flags %>% filter(state %in% c("queue_fresh", "queue_retry"))
  if (nrow(q) == 0) {
    message("nothing to queue")
    return(invisible(flags))
  }

  stopifnot(!any(is.na(q[["p1_name"]])), !any(is.na(q[["p2_name"]])),
            all(nzchar(q[["p1_name"]])), all(nzchar(q[["p2_name"]])),
            all(fs::dir_exists(q[["afpd_dir"]])),
            all(q[["n_iptm"]] == E), all(q[["num_ark"]] == E),
            identical(unique(as.character(fs::path_file(q[["run_dir"]]))),
                      tree[["run_name"]]))

  ## A second pass while the first array is still writing metrics_v2.csv into
  ## these same directories would have two workers writing the same files.
  inflight <- squeue_extract_metrics()
  if (length(inflight)) {
    msg <- paste0(length(inflight), " extract_metrics job(s) still in flight: ",
                  paste(inflight, collapse = ", "))
    if (submit) stop("refusing to submit -- ", msg, call. = FALSE)
    message("!! ", msg, " -- the state above may be mid-write")
  }

  job_dir <- fs::path(JOB_ROOT, paste0(STAMP, "_", key))
  submit_metric_jobs(q %>% select(afpd_dir, p1_name, p2_name, run_dir),
                     group_size = GROUP_SIZE,
                     job_dir    = job_dir,
                     dry_run    = !submit)
  invisible(flags)
}

squeue_extract_metrics <- function() {
  user <- Sys.getenv("USER", Sys.info()[["user"]])
  out <- tryCatch(
    suppressWarnings(
      system2("squeue", c("-h", "-u", user, "-n", "extract_metrics", "-o", "%i"),
              stdout = TRUE)),
    error = function(e) structure(character(0), status = -1L,
                                  reason = conditionMessage(e)))
  st <- attr(out, "status")
  if (!is.null(st) && st != 0) {
    stop("could not run squeue to check for in-flight extract_metrics jobs (",
         paste(c(paste("exit", st), attr(out, "reason")), collapse = ": "),
         ") -- check by hand", call. = FALSE)
  }
  out[nzchar(out)]
}

## --- dispatch -------------------------------------------------------------

main <- function(args) {
  usage <- paste0("usage: Rscript afpd_prepare_metric_jobs.R <",
                  paste(names(TREES), collapse = "|"),
                  "> <survey|rename|queue> [--submit]")
  if (length(args) < 2 || !args[1] %in% names(TREES) ||
      !args[2] %in% c("survey", "rename", "queue")) {
    stop(usage, call. = FALSE)
  }
  key    <- args[1]
  step   <- args[2]
  submit <- "--submit" %in% args
  if (submit && step != "queue") stop("--submit only applies to 'queue'", call. = FALSE)
  tree <- TREES[[key]]
  stopifnot(dir.exists(tree[["run_dir"]]),
            identical(as.character(fs::path_file(tree[["run_id_dir"]])),
                      tree[["run_name"]]))
  switch(step,
         survey = step_survey(tree, key),
         rename = step_rename(tree, key),
         queue  = step_queue(tree, key, submit = submit))
  invisible(NULL)
}

if (sys.nframe() == 0L) main(commandArgs(trailingOnly = TRUE))
