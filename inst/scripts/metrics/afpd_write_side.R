## ---------------------------------------------------------------------------
## afpd_write_side.R  --  raw AlphaPulldown directories -> metric jobs queued
##
## Three trees, one parameter table, and the six rename stages driven
## individually so an assertion sits in front of every irreversible action and
## a failure cannot leave a half-renamed tree behind.
##
## WHERE TO RUN IT.  Not on a login node.  lf_gate() calls afpd_check_metrics(),
## which materialises every metrics_v2.csv into a list column at once, and
## lf_stage5() reads and rewrites the 43M-row random_codes.csv.  Both get
## OOM-killed on a login node.  Ask for:
##
##   salloc -p normal -c 8 --mem=64G -t 4:00:00
##   ml R/4.1.2
##   R
##
## THE SEQUENCE, one step at a time, read between each:
##
##   source(system.file("scripts/metrics/afpd_write_side.R", package = "fap"))
##   lf_guard()
##   tr <- LF_TREES$A                              # then $B, then $C
##   ce <- lf_census(tr);            print(ce %>% count(kind, status))
##   qu <- lf_quarantine(tr, ce)                   # IRREVERSIBLE (moves dirs)
##   s1 <- lf_stage1_parse_dirnames(tr, ce)        # pure
##   s2 <- lf_stage2_make_new_dirnames(tr, s1)     # pure
##   s3 <- lf_stage3_rename_dirs(tr, s2)           # IRREVERSIBLE
##   s4 <- lf_stage4_parse_files(tr, s3)
##   s5 <- lf_stage5_make_file_names(tr, s4)       # consumes random codes
##   s6 <- lf_stage6_rename_files(tr, s5)          # IRREVERSIBLE
##   g  <- lf_gate(tr);              print(g$summary)
##   lf_submit(tr, g)                              # dry run: writes group files
##   lf_submit(tr, g, dry_run = FALSE)             # IRREVERSIBLE: sbatch
##
## The script is installed with the package, so either of these gets it:
##   source(system.file("scripts/metrics/afpd_write_side.R", package = "fap"))
##   source("/home/groups/ebutcher/kevin/ligandFinder/afpd_write_side.R")
##
## Run the same sequence again, verbatim, after the --mem=64G restart lands.
## Every step is a no-op on work already done; only genuinely new directories
## move.  Nothing is ever marked finished on the strength of a short
## prediction set.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr); library(purrr); library(tibble); library(tidyr)
  library(stringr); library(fs); library(data.table); library(jsonlite)
})
if (!requireNamespace("ligandFinder", quietly = TRUE))
  stop("ligandFinder is not installed on this libPath")
library(fap)

LF_RC_FILE  <- "/home/groups/ebutcher/kevin/ligandFinder/random_codes.csv"
LF_JOB_ROOT <- "/scratch/groups/ebutcher/deorphan/ligandFinder/metric_jobs"

## A renamed directory is hRECEPTOR_hLIGANDxSTARTxEND and must hold exactly one
## "_": afpd_check_metrics() -> parse_proteins() derives num_proteins from
## max(str_count(name, "_")) across the batch, so one two-underscore name would
## silently re-parse every other directory in the tree.
LF_RENAMED_RE <- "^h[A-Z0-9]+_h[A-Z0-9]+x[0-9]+x[0-9]+$"
LF_RAW_RE     <- "_and_"

lf_tree <- function(label, dir, n_pred, relax, n_intended) {
  n_xtr <- as.integer(n_pred) * (1L + as.integer(relax))
  list(label = label, dir = dir,
       n_pred     = as.integer(n_pred),   # predictions == ark == con == pae == length(iptm)
       relax      = relax,
       n_xtr      = n_xtr,                # unrelaxed (+ relaxed) pdbs
       n_models   = max(as.integer(n_pred), n_xtr),  # afpd_check_metrics' max(ark, xtr)
       n_intended = as.integer(n_intended),
       ## submit_metric_jobs() only ever takes fs::path_file(run_dir), to build
       ## the run_id field of the group files.  Run B's directory is
       ## .../cxl14_75Y/models, whose basename is "models", so the label travels
       ## as a bare string instead of a path.
       work       = fs::path(dir, "_lf_pipeline"),
       quarantine = fs::path(dir, "_quarantine"))
}

LF_TREES <- list(
  A = lf_tree("cxl14_x75", "/scratch/groups/ebutcher/deorphan/models/cxl14_x75",
              n_pred = 5,  relax = FALSE, n_intended = 288),
  B = lf_tree("cxl14_75Y", "/oak/stanford/groups/ebutcher/kevin/cxl14_75Y/models",
              n_pred = 25, relax = TRUE,  n_intended = 4),
  C = lf_tree("vtm2a_157", "/scratch/groups/ebutcher/deorphan/models/vtm2a_157",
              n_pred = 5,  relax = FALSE, n_intended = 576)
)

## --- helpers ---------------------------------------------------------------

lf_hr <- function(...) message("\n== ", ..., " ", strrep("=", 24))

lf_assert <- function(ok, msg, offenders = NULL) {
  ok <- ok %in% TRUE                     # NA counts as failure, never as pass
  if (all(ok)) return(invisible(TRUE))
  if (!is.null(offenders))
    message("offenders (up to 20 of ", sum(!ok), "): ",
            paste(utils::head(unique(offenders[!ok]), 20), collapse = ", "))
  stop("ASSERTION FAILED: ", msg, call. = FALSE)
}

lf_guard <- function() {
  if (!nzchar(Sys.getenv("SLURM_JOB_ID")))
    stop("No SLURM_JOB_ID.\n",
         "lf_gate() holds every metrics_v2.csv in memory at once and lf_stage5()\n",
         "rewrites a 43M-row CSV; both are OOM-killed on a login node. Run:\n",
         "  salloc -p normal -c 8 --mem=64G -t 4:00:00")
  cpus <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
  if (is.na(cpus) || cpus < 1L) cpus <- 1L
  workers <- max(1L, cpus %/% 2L)
  future::plan(strategy = future::multicore(workers = workers))
  options(future.globals.maxSize = 32 * 1024^3)
  message("job ", Sys.getenv("SLURM_JOB_ID"), ": ", cpus, " cpus, ", workers,
          " future workers, future.globals.maxSize 32 GiB")
  invisible(workers)
}

## gpcr_list comes from R/sysdata.rda and namespace lookup wins over any
## readRDS(system.file("extdata/gpcr_list.rds")).  Four of the 845 rows hold a
## bare list rather than a data frame, so nrow() returns NULL there; 455 more
## hold a zero-row table.  459 receptors can therefore never produce depth.
lf_bw_rows <- function() {
  gl <- fap:::gpcr_list
  n <- vapply(gl[["bw: full_table"]],
              function(z) { k <- nrow(z); if (is.null(k)) NA_integer_ else as.integer(k) },
              integer(1))
  tibble(p1_name = gl[["uniprot_name"]], bw_rows = n) %>%
    distinct(p1_name, .keep_all = TRUE)
}

## Count unused random codes without pulling 43M rows into R.
lf_unused_codes <- function(rc_file = LF_RC_FILE) {
  lf_assert(file.exists(rc_file), paste("random codes file missing:", rc_file))
  prog <- tempfile(fileext = ".awk")
  on.exit(unlink(prog), add = TRUE)
  writeLines('NR>1 { c[$NF]++ } END { for (k in c) printf "%s\\t%d\\n", k, c[k] }', prog)
  out <- suppressWarnings(system2("awk", c("-F,", "-f", shQuote(prog), shQuote(rc_file)),
                                  stdout = TRUE, stderr = FALSE))
  n <- NA_integer_
  if (length(out)) {
    kv <- do.call(rbind, strsplit(out, "\t", fixed = TRUE))
    if (!is.null(kv) && ncol(kv) == 2L && "unused" %in% kv[, 1])
      n <- as.integer(kv[kv[, 1] == "unused", 2][1])
  }
  if (is.na(n)) {                       # no awk: read the one column
    u <- data.table::fread(rc_file, select = "usage", showProgress = FALSE)
    n <- sum(u[["usage"]] == "unused")
  }
  message("unused random codes: ", format(n, big.mark = ","))
  as.integer(n)
}

## ranking_debug.json, read with the predicates the pipeline actually imposes:
## parse_afpd_files() runs as_tibble(fromJSON(.)) %>% select(order), which needs
## "order" present and all three fields the same length; import_raw_metrics()
## needs names(iptm) and iptm+ptm.  A file that fails any of these aborts
## stage 4 -- after stage 3 has already renamed the whole tree.
lf_read_rd <- function(d) {
  out <- tibble(rd_present = FALSE, rd_valid = FALSE, rd_iptm_named = FALSE,
                n_order = NA_integer_, n_iptm = NA_integer_, n_iptmptm = NA_integer_)
  f <- file.path(d, "ranking_debug.json")
  if (!file.exists(f)) return(out)
  out$rd_present <- TRUE
  j <- tryCatch(jsonlite::read_json(f), error = function(e) NULL)
  if (is.null(j) || !is.list(j)) return(out)
  nm <- names(j)
  len <- function(k) if (k %in% nm) length(j[[k]]) else NA_integer_
  out$n_order   <- len("order")
  out$n_iptm    <- len("iptm")
  out$n_iptmptm <- len("iptm+ptm")
  out$rd_iptm_named <- "iptm" %in% nm && !is.null(names(j[["iptm"]])) &&
                       all(nzchar(names(j[["iptm"]])))
  out$rd_valid <- all(c("order", "iptm", "iptm+ptm") %in% nm) &&
                  !any(is.na(c(out$n_order, out$n_iptm, out$n_iptmptm))) &&
                  out$n_order == out$n_iptm && out$n_iptm == out$n_iptmptm &&
                  out$n_order > 0L && out$rd_iptm_named
  out
}

lf_file_counts <- function(d) {
  f <- list.files(d)
  tibble(
    raw_ranked    = sum(grepl("^ranked_[0-9]+\\.pdb$", f)),
    raw_con       = sum(grepl("^confidence_model_.*\\.json$", f)),
    raw_pae       = sum(grepl("^pae_model_.*\\.json$", f)),
    raw_unrelaxed = sum(grepl("^unrelaxed_model_.*\\.pdb$", f)),
    raw_relaxed   = sum(grepl("^relaxed_model_.*\\.pdb$", f)),
    ## the regexes import_raw_metrics() and afpd_check_metrics() grep with
    n_ark = sum(grepl("_ark_.*\\.pdb$", f)),
    n_con = sum(grepl("_con_.*\\.json$", f)),
    n_pae = sum(grepl("_pae_.*\\.json$", f)),
    n_xtr = sum(grepl("_xtr_.*\\.pdb$", f)),
    has_log = "file_name_log.csv" %in% f)
}

lf_write_tsv <- function(x, path) {
  fs::dir_create(fs::path_dir(path))
  data.table::fwrite(x %>% select(!where(is.list)), file = path, sep = "\t")
  message("wrote ", path, " (", nrow(x), " rows)")
  invisible(path)
}

lf_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

## === STEP 1  census (read-only) ===========================================

lf_census <- function(tree) {
  lf_hr("census ", tree$label)
  lf_assert(fs::dir_exists(tree$dir), paste("run tree does not exist:", tree$dir))
  ## The rlx token in a renamed file name is a single "_u_" / "_r_", and
  ## import_raw_metrics() finds it with str_extract(.., "_[ru]_") -- first match
  ## wins, so the run label must not contain one.
  lf_assert(!grepl("_[ru]_", paste0("_", tree$label, "_")),
            paste0("run label '", tree$label, "' contains _r_ or _u_, which would ",
                   "shadow the rlx token"))

  all_dirs <- fs::dir_ls(tree$dir, type = "directory", recurse = FALSE)
  nm <- fs::path_file(all_dirs)
  ctrl <- startsWith(nm, "_") | startsWith(nm, ".")
  if (any(ctrl)) message(sum(ctrl), " control directories skipped: ",
                         paste(nm[ctrl], collapse = ", "))

  ce <- tibble(tree_label = tree$label,
               dir_name = nm[!ctrl],
               afpd_dir = as.character(all_dirs[!ctrl])) %>%
    mutate(kind = case_when(grepl(LF_RAW_RE, dir_name)     ~ "raw",
                            grepl(LF_RENAMED_RE, dir_name) ~ "renamed",
                            TRUE                           ~ "unknown"))

  ## An empty tree is normal: a fresh run directory, or one whose every entry
  ## was quarantined. bind_rows() over nothing yields a tibble with no columns,
  ## so take the schema from a probe of a directory that holds none of these
  ## files.
  if (nrow(ce) == 0) {
    ce <- ce %>% bind_cols(lf_read_rd(tempdir())[0, ], lf_file_counts(tempdir())[0, ]) %>%
            mutate(raw_full = logical(0), status = character(0))
    message("no pairing directories in this tree; every step below is a no-op")
    return(invisible(ce))
  }

  ce <- bind_cols(ce,
                  bind_rows(furrr::future_map(ce$afpd_dir, lf_read_rd)),
                  bind_rows(furrr::future_map(ce$afpd_dir, lf_file_counts)))

  raw_full <- with(ce, raw_ranked == tree$n_pred & raw_con == tree$n_pred &
                        raw_pae == tree$n_pred & raw_unrelaxed == tree$n_pred &
                        (!tree$relax | raw_relaxed == tree$n_pred))

  ce <- ce %>% mutate(
    raw_full = raw_full,
    status = case_when(
      kind == "unknown"                                 ~ "unknown_name",
      kind == "renamed" & n_ark == 0 & raw_ranked > 0    ~ "half_renamed",
      kind == "renamed"                                 ~ "renamed",
      !rd_present                                       ~ "pending",
      !rd_valid                                         ~ "malformed",
      n_iptm != tree$n_pred | !raw_full                  ~ "short",
      TRUE                                              ~ "ready"))

  print(ce %>% count(kind, status))
  message("intended pairs: ", tree$n_intended, " | directories present: ", nrow(ce))
  if (any(ce$status == "unknown_name"))
    message("unknown_name: left completely alone, never renamed and never queued.")
  if (any(ce$status == "pending"))
    message("pending: no ranking_debug.json yet (still running, or died before ranking).\n",
            "Left raw so the restart can overwrite them in place.")
  if (any(ce$status == "short"))
    message("short: fewer than ", tree$n_pred, " predictions. Left RAW on purpose -- ",
            "renaming them\nwould make the restart's fresh output collide with the renamed copy.")
  if (any(ce$status == "half_renamed"))
    message("half_renamed: directory renamed but files not. Repair with ",
            "lf_resume_files(tree, ce)\nbefore lf_gate().")
  invisible(ce)
}

## === STEP 2  quarantine malformed ranking_debug.json (IRREVERSIBLE) =======

lf_quarantine <- function(tree, ce) {
  lf_hr("quarantine ", tree$label)
  lf_assert(all(ce$tree_label == tree$label), "census is from a different tree")

  bad <- ce %>% filter(status == "malformed")
  if (nrow(bad) == 0) { message("nothing malformed; no-op"); return(invisible(bad)) }

  from <- bad$afpd_dir
  to   <- fs::path(tree$quarantine, bad$dir_name)

  lf_assert(all(grepl(LF_RAW_RE, bad$dir_name)), "would quarantine a non-raw directory", bad$dir_name)
  lf_assert(all(bad$rd_present & !bad$rd_valid), "candidate is not actually malformed", bad$dir_name)
  lf_assert(all(fs::dir_exists(from)), "source directory vanished", from)
  lf_assert(!any(fs::file_exists(to)), "quarantine target already exists", to)

  fs::dir_create(tree$quarantine)
  lf_write_tsv(bad, fs::path(tree$work, paste0("quarantined_", lf_stamp(), ".tsv")))

  ok <- file.rename(from, to)
  lf_assert(ok, "file.rename() of a quarantine candidate returned FALSE", from)
  lf_assert(all(fs::dir_exists(to)) && !any(fs::dir_exists(from)), "quarantine move did not land")
  message("quarantined ", nrow(bad), " directories to ", tree$quarantine,
          "\n(the name starts with '_' and holds no '_and_', so neither lf_census() nor\n",
          "parse_dirname(afpd_raw = TRUE) will ever see it again)")
  invisible(bad)
}

## === STAGE 1  parse_dirname (pure: no filesystem effect) ==================

lf_stage1_parse_dirnames <- function(tree, ce) {
  lf_hr("stage 1 parse_dirname ", tree$label)
  lf_assert(all(ce$tree_label == tree$label), "census is from a different tree")

  ready <- ce %>% filter(status == "ready")
  if (nrow(ready) == 0) {
    message("no 'ready' raw directories; the rename half of this pass is a no-op")
    return(invisible(NULL))
  }
  message("parsing ", nrow(ready), " names")

  ## pairing_dir keeps this pure. The delimiter defaults ("_and_" / "_" / "-")
  ## are the ones do_renaming() uses and the ones the measured raw shape needs.
  pd <- parse_dirname(pairing_dir = ready$dir_name)

  lf_assert(nrow(pd) == nrow(ready), "parse_dirname changed the row count")
  lf_assert(setequal(pd$afpd_dir_name, ready$dir_name), "parse_dirname changed the name set")

  ## An accession the Entry -> Entry Name converter cannot resolve becomes NA,
  ## and filter(value != "") then drops it silently, so model_name is built from
  ## fewer parts instead of failing.
  need <- c("p1.name", "p1.id", "p2.name", "p2.id", "p2.range")
  have <- map(pd$parsed_pair, \(x) paste(x$protein, x$annotation, sep = "."))
  lf_assert(map_lgl(have, \(v) all(need %in% v)),
            "a directory has an unresolvable accession or no ligand range", pd$afpd_dir_name)
  lf_assert(!is.na(pd$model_name) & nzchar(pd$model_name), "empty model_name", pd$afpd_dir_name)
  lf_assert(!is.na(pd$model_id)   & nzchar(pd$model_id),   "empty model_id",   pd$afpd_dir_name)

  pick <- function(x, p, a) {
    v <- x$value[x$protein == p & x$annotation == a]
    if (length(v)) v[1] else NA_character_
  }
  pd <- pd %>% mutate(
    p1_name = map_chr(parsed_pair, pick, "p1", "name"),
    p2_name = map_chr(parsed_pair, pick, "p2", "name"),
    ## make_model_names() appends gpcr_list's first_AA-last_AA, and
    ## make_new_dirname(exclude_p1_range = TRUE) then drops the SUBMITTED p1
    ## range from the directory name. import_raw_metrics() re-parses the RENAMED
    ## name, so metrics_v2.csv will carry gpcr_list's range and not what was
    ## actually run. These two columns are the only surviving record of it.
    submitted_p1_range = map_chr(parsed_pair, pick, "p1", "range"),
    submitted_p2_range = map_chr(parsed_pair, pick, "p2", "range")) %>%
    left_join(ready %>% select(afpd_dir_name = dir_name, afpd_dir, n_iptm),
              by = "afpd_dir_name")

  print(pd %>% count(p1_name_has_submitted_range = !is.na(submitted_p1_range) &
                       nzchar(submitted_p1_range), p2_name, submitted_p2_range))
  invisible(pd)
}

## === STAGE 2  make_new_dirname (pure) =====================================

lf_stage2_make_new_dirnames <- function(tree, s1) {
  lf_hr("stage 2 make_new_dirname ", tree$label)
  if (is.null(s1)) { message("nothing to plan"); return(invisible(NULL)) }

  nd <- make_new_dirname(input = s1,
                         delim_proteins = "_", delim_ranges = "x", delim_start_end = "x",
                         p1_prefix = "h", p1_suffix = NA, p2_prefix = "h",
                         exclude_p1_range = TRUE)
  ## make_new_dirname() keeps only its own columns; put the provenance back.
  nd <- nd %>% left_join(s1 %>% select(afpd_dir_name, afpd_dir, p1_name, p2_name,
                                       submitted_p1_range, submitted_p2_range, n_iptm),
                         by = "afpd_dir_name")

  lf_assert(nrow(nd) == nrow(s1), "make_new_dirname changed the row count")
  lf_assert(!is.na(nd$new_dir_name) & nzchar(nd$new_dir_name), "empty new_dir_name", nd$afpd_dir_name)
  lf_assert(grepl(LF_RENAMED_RE, nd$new_dir_name),
            paste0("new_dir_name does not match ", LF_RENAMED_RE,
                   "; afpd_check_metrics() would mis-parse the tree"), nd$new_dir_name)
  lf_assert(str_count(nd$new_dir_name, "_") == 1L,
            "new_dir_name has other than one underscore", nd$new_dir_name)
  lf_assert(nd$new_dir_name != nd$afpd_dir_name, "new name equals old name", nd$afpd_dir_name)

  ## The collision exclude_p1_range = TRUE creates: a bare receptor and a
  ## range-truncated receptor against the same ligand window map onto one name.
  dups <- nd %>% count(new_dir_name) %>% filter(n > 1)
  if (nrow(dups))
    print(nd %>% filter(new_dir_name %in% dups$new_dir_name) %>%
            select(afpd_dir_name, new_dir_name, model_name), n = Inf)
  lf_assert(nrow(dups) == 0,
            "two raw directories collide on one new_dir_name (bare vs range-truncated receptor)",
            dups$new_dir_name)

  on_disk <- fs::path_file(fs::dir_ls(tree$dir, type = "directory", recurse = FALSE))
  clash <- intersect(nd$new_dir_name, on_disk)
  lf_assert(length(clash) == 0,
            paste("a planned new_dir_name already exists on disk. file.rename() onto an",
                  "EMPTY existing directory SUCCEEDS and replaces it, so this is a",
                  "data-loss path, not a harmless no-op", sep = "\n"), clash)

  ## A receptor absent from gpcr_list, or holding an empty 'bw: full_table',
  ## can never produce depth. Renaming is still correct; lf_gate() excludes it.
  nobw <- nd %>% left_join(lf_bw_rows(), by = "p1_name") %>%
            filter(is.na(bw_rows) | bw_rows == 0)
  if (nrow(nobw))
    message("NOTE: ", nrow(nobw), " of ", nrow(nd), " receptors have no usable ",
            "'bw: full_table'.\nlf_gate() will report them as incapable_no_bw and never queue them.")

  nd <- nd %>% mutate(tree_label = tree$label, .before = 1)
  path <- fs::path(tree$work, paste0("rename_plan_", lf_stamp(), ".tsv"))
  lf_write_tsv(nd %>% select(tree_label, afpd_dir_name, new_dir_name, model_name, model_id,
                             p1_name, p2_name, submitted_p1_range, submitted_p2_range, n_iptm),
               path)
  attr(nd, "plan_path") <- path
  print(nd %>% select(afpd_dir_name, new_dir_name) %>% head(5))
  invisible(nd)
}

## === STAGE 3  rename_dir (IRREVERSIBLE) ===================================

lf_stage3_rename_dirs <- function(tree, s2) {
  lf_hr("stage 3 rename_dir ", tree$label, "  IRREVERSIBLE")
  if (is.null(s2)) { message("nothing to rename"); return(invisible(NULL)) }

  from <- fs::path(tree$dir, s2$afpd_dir_name)
  to   <- fs::path(tree$dir, s2$new_dir_name)

  lf_assert(!is.null(attr(s2, "plan_path")) && fs::file_exists(attr(s2, "plan_path")),
            "the stage 2 plan TSV is not on disk; without it a crash mid-rename is unrecoverable")
  lf_assert(fs::dir_exists(from), "source directory vanished since stage 2", from)
  lf_assert(!fs::file_exists(to), "target name appeared since stage 2", to)
  lf_assert(anyDuplicated(s2$new_dir_name) == 0, "duplicate target names")

  ## stage 5 calls get_codes(), which stop()s when codes run short. That is
  ## AFTER this rename, so prove the codes exist before touching the tree.
  unused <- lf_unused_codes()
  lf_assert(unused >= nrow(s2),
            paste0("only ", unused, " unused random codes for ", nrow(s2),
                   " directories; stage 5 would stop() after this rename"))

  res <- rename_dir(run_dir = tree$dir, input = s2,
                    from = "afpd_dir_name", to = "new_dir_name")

  lf_assert(res$rename_status, "file.rename() returned FALSE for a directory", s2$afpd_dir_name)
  lf_assert(fs::dir_exists(to),   "target directory missing after rename", to)
  lf_assert(!fs::dir_exists(from), "source directory still present after rename", from)
  message(nrow(s2), " directories renamed")
  invisible(res)
}

## === STAGE 4  parse_afpd_files ============================================

lf_stage4_parse_files <- function(tree, s3) {
  lf_hr("stage 4 parse_afpd_files ", tree$label)
  if (is.null(s3)) { message("nothing to parse"); return(invisible(NULL)) }

  pf <- parse_afpd_files(input = s3, dir_name = "new_dir_name", run_dir = tree$dir)

  ## parse_afpd_files() silently filter()s out directories whose JSON will not
  ## parse. Those are now renamed but unprocessed, so name them loudly.
  lost <- setdiff(s3$new_dir_name, pf$new_dir_name)
  lf_assert(length(lost) == 0,
            paste("parse_afpd_files() dropped a directory: its ranking_debug.json went bad",
                  "after the census. It is renamed but unprocessed -- repair with",
                  "lf_resume_files(tree, lf_census(tree))", sep = "\n"), lost)

  lf_assert(map_int(pf$ranks, nrow) == tree$n_pred,
            paste0("the ranks table is not ", tree$n_pred, " rows"), pf$new_dir_name)

  cnt <- pf %>% mutate(
    n_ranked    = map_int(files, \(x) sum(x$file_type == "ranked", na.rm = TRUE)),
    n_con       = map_int(files, \(x) sum(x$file_type == "confidence", na.rm = TRUE)),
    n_pae       = map_int(files, \(x) sum(x$file_type == "pae", na.rm = TRUE)),
    n_unrelaxed = map_int(files, \(x) sum(x$file_type == "unrelaxed", na.rm = TRUE)),
    n_relaxed   = map_int(files, \(x) sum(x$file_type == "relaxed", na.rm = TRUE)),
    n_badtype   = map_int(files, \(x) {
      ft <- x$file_type[!is.na(x$model)]
      sum(!ft %in% names(fap:::file_type_conv)) }),
    n_norank    = map_int(files, \(x) sum(!is.na(x$model) & is.na(x$rank))))

  lf_assert(cnt$n_ranked    == tree$n_pred, "ranked file count wrong", cnt$new_dir_name)
  lf_assert(cnt$n_con       == tree$n_pred, "confidence file count wrong", cnt$new_dir_name)
  lf_assert(cnt$n_pae       == tree$n_pred, "pae file count wrong", cnt$new_dir_name)
  lf_assert(cnt$n_unrelaxed == tree$n_pred, "unrelaxed file count wrong", cnt$new_dir_name)
  lf_assert(cnt$n_relaxed == (if (tree$relax) tree$n_pred else 0L),
            paste0("relaxed file count is not ", if (tree$relax) tree$n_pred else 0L,
                   "; this tree's driver ", if (tree$relax) "does" else "does not",
                   " pass --models_to_relax"), cnt$new_dir_name)
  lf_assert(cnt$n_badtype == 0L,
            "a model-bearing file has a file_type absent from file_type_conv; its new name\nwould silently lose the type token",
            cnt$new_dir_name)
  lf_assert(cnt$n_norank == 0L,
            "a model-bearing file did not join to ranking_debug.json's order", cnt$new_dir_name)

  message("parsed ", nrow(pf), " directories, ",
          sum(map_int(pf$files, \(x) sum(!is.na(x$model)))), " model-bearing files")
  invisible(pf)
}

## === STAGE 5  make_new_file_names (consumes random codes) =================

lf_stage5_make_file_names <- function(tree, s4) {
  lf_hr("stage 5 make_new_file_names ", tree$label)
  if (is.null(s4)) { message("nothing to name"); return(invisible(NULL)) }

  unused <- lf_unused_codes()
  lf_assert(unused >= nrow(s4),
            paste0("only ", unused, " unused codes for ", nrow(s4), " directories"))
  message("reading and rewriting ", LF_RC_FILE, " (43M rows); allow a few minutes")

  mf <- make_new_file_names(input = s4, dir_name = "new_dir_name",
                            run_name = tree$label, site = "SU", submitter = "KB",
                            algorithm = "AF2v3", random_seed = 42, rc_file = LF_RC_FILE)

  ## import_raw_metrics() extracts "_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+".
  ## final_code is site + submitter + random_code + model_code = 4 upper then
  ## 7 lower, which only holds while get_codes() is inside the letters-only
  ## block of random_codes.csv.
  lf_assert(grepl("^[a-z]{5}$", mf$random_code),
            "get_codes() returned a code outside the letters-only block; the file names\nwould not match _[A-Z]{4}[a-z]{7}_",
            mf$random_code)
  lf_assert(anyDuplicated(mf$random_code) == 0, "duplicate random code in this batch")

  ## The rlx token is present on ark/xtr names and absent on con/pae names, so
  ## the separator before the code is optional.
  code <- "[A-Z]{4}[a-z]{7}_s[0-9]+m[0-9]+p[0-9]+"
  ## The two extractions import_raw_metrics() actually performs on an ark name.
  irm_code <- "_[A-Z]{4}[a-z]{7}_s[0-9]+m[0-9]+p[0-9]+_r[0-9]+"
  irm_rlx  <- "_[ru]_"
  chk <- mf %>% mutate(
    n_na  = map_int(files, \(x) sum(is.na(x$new_file_name))),
    n_dup = map_int(files, \(x) sum(duplicated(x$new_file_name))),
    ark   = map_int(files, \(x) sum(grepl(paste0("_ark_(.+_)?", code, "_r[0-9]+\\.pdb$"), x$new_file_name))),
    con   = map_int(files, \(x) sum(grepl(paste0("_con_(.+_)?", code, "(_r[0-9]+)?\\.json$"), x$new_file_name))),
    pae   = map_int(files, \(x) sum(grepl(paste0("_pae_(.+_)?", code, "(_r[0-9]+)?\\.json$"), x$new_file_name))),
    xtr   = map_int(files, \(x) sum(grepl(paste0("_xtr_(.+_)?", code, "(_r[0-9]+)?\\.pdb$"), x$new_file_name))),
    irm   = map_int(files, \(x) {
      a <- grep("_ark_.*\\.pdb$", x$new_file_name, value = TRUE)
      sum(grepl(irm_code, a) & grepl(irm_rlx, a)) }),
    clash = map2_int(new_dir_name, files, \(d, x) {
      have <- list.files(fs::path(tree$dir, d))
      sum(x$new_file_name[x$new_file_name != x$og_file_name] %in% have) }))

  lf_assert(chk$n_na  == 0L, "NA new_file_name", chk$new_dir_name)
  lf_assert(chk$n_dup == 0L, "two files would be renamed onto one name", chk$new_dir_name)
  lf_assert(chk$ark == tree$n_pred, "ark names do not match import_raw_metrics' pattern", chk$new_dir_name)
  lf_assert(chk$con == tree$n_pred, "con names do not match the _con_ pattern", chk$new_dir_name)
  lf_assert(chk$pae == tree$n_pred, "pae names do not match the _pae_ pattern", chk$new_dir_name)
  lf_assert(chk$xtr == tree$n_xtr,  "xtr names do not match the _xtr_ pattern", chk$new_dir_name)
  lf_assert(chk$irm == tree$n_pred,
            "an ark name fails the exact extractions import_raw_metrics() performs:\nstr_extract(., \"_[A-Z]{4}[a-z]{7}_s\\\\d+m\\\\d+p\\\\d+_r\\\\d+\") and str_extract(., \"_[ru]_\")",
            chk$new_dir_name)
  lf_assert(chk$clash == 0L, "a planned file name already exists in the directory", chk$new_dir_name)

  st <- lf_stamp()
  path <- fs::path(tree$work, paste0("file_plan_", st, ".tsv"))
  lf_write_tsv(mf %>% mutate(tree_label = tree$label) %>%
                 select(tree_label, new_dir_name, random_code), path)
  lf_write_tsv(mf %>% select(new_dir_name, files) %>% tidyr::unnest(files) %>%
                 select(new_dir_name, og_file_name, new_file_name, file_type, rlx,
                        model_num, pred_num, rank),
               fs::path(tree$work, paste0("file_plan_detail_", st, ".tsv")))
  attr(mf, "plan_path") <- path
  invisible(mf)
}

## === STAGE 6  rename_files (IRREVERSIBLE) =================================

lf_stage6_rename_files <- function(tree, s5) {
  lf_hr("stage 6 rename_files ", tree$label, "  IRREVERSIBLE")
  if (is.null(s5)) { message("nothing to rename"); return(invisible(NULL)) }

  lf_assert(!is.null(attr(s5, "plan_path")) && fs::file_exists(attr(s5, "plan_path")),
            "the stage 5 plan TSV is not on disk")
  lf_assert(fs::dir_exists(fs::path(tree$dir, s5$new_dir_name)),
            "a target directory vanished since stage 5", s5$new_dir_name)
  ## file_name_log.csv means this directory already went through stage 6;
  ## re-running would rename nothing and overwrite the log.
  already <- fs::file_exists(fs::path(tree$dir, s5$new_dir_name, "file_name_log.csv"))
  lf_assert(!already, "directory already has file_name_log.csv (stage 6 already ran here)",
            s5$new_dir_name)

  miss <- map2_int(s5$new_dir_name, s5$files, \(d, x)
    sum(!fs::file_exists(fs::path(tree$dir, d, x$og_file_name))))
  lf_assert(miss == 0L, "a source file named in the stage 5 plan is missing", s5$new_dir_name)

  rename_files(run_dir = tree$dir, input = s5, dir_name = "new_dir_name",
               from = "og_file_name", to = "new_file_name")

  post <- bind_rows(map(fs::path(tree$dir, s5$new_dir_name), lf_file_counts)) %>%
            mutate(new_dir_name = s5$new_dir_name)
  lf_assert(post$n_ark == tree$n_pred, "post-rename ark count wrong", post$new_dir_name)
  lf_assert(post$n_con == tree$n_pred, "post-rename con count wrong", post$new_dir_name)
  lf_assert(post$n_pae == tree$n_pred, "post-rename pae count wrong", post$new_dir_name)
  lf_assert(post$n_xtr == tree$n_xtr,  "post-rename xtr count wrong", post$new_dir_name)
  lf_assert(post$raw_ranked == 0L, "raw ranked_*.pdb still present after rename", post$new_dir_name)
  lf_assert(post$has_log, "file_name_log.csv was not written", post$new_dir_name)
  message(nrow(post), " directories fully renamed")
  invisible(post)
}

## --- repair path: directory renamed, files not ----------------------------

lf_resume_files <- function(tree, ce) {
  lf_hr("resume stages 4-6 ", tree$label)
  hr <- ce %>% filter(status == "half_renamed")
  if (nrow(hr) == 0) { message("nothing half-renamed"); return(invisible(NULL)) }
  lf_assert(hr$rd_valid, "half-renamed directory has a malformed ranking_debug.json", hr$dir_name)
  s4 <- lf_stage4_parse_files(tree, tibble(new_dir_name = hr$dir_name))
  s5 <- lf_stage5_make_file_names(tree, s4)
  lf_stage6_rename_files(tree, s5)
}

## === STEP 3  the queue gate (HEAVY -- needs the allocation) ===============

lf_gate <- function(tree, max_attempts = 2L) {
  lf_hr("gate ", tree$label)
  lf_assert(nzchar(Sys.getenv("SLURM_JOB_ID")),
            "lf_gate() must not run on a login node; call lf_guard() first")

  nm <- fs::path_file(fs::dir_ls(tree$dir, type = "directory", recurse = FALSE))
  nm <- nm[grepl(LF_RENAMED_RE, nm)]
  lf_assert(length(nm) > 0, "no renamed directories in this tree")
  lf_assert(str_count(nm, "_") == 1L, "a renamed directory name has other than one underscore", nm)

  runs <- tibble(afpd_dir_name = nm,
                 afpd_dir = as.character(fs::path(tree$dir, nm)),
                 run_dir = tree$label)   # submit_metric_jobs() takes path_file() of this

  message("afpd_check_metrics() over ", nrow(runs), " directories; this is the memory peak")
  cm <- afpd_check_metrics(runs)

  ## Reduce the metrics list column to flags and drop it in the same pipeline.
  ## metrics_good from afpd_check_metrics() is not used: it tests
  ## nrow(metrics) == num_models, and for run B num_models is max(ark, xtr) = 50
  ## while only 25 predictions exist. afpd_check_contacts() is deliberately not
  ## called at all -- it readRDS()s every metrics_v2c.rds, multi-GB over this
  ## many directories, and tells us nothing about whether a job is needed.
  flags <- cm %>%
    mutate(metrics_rows      = map_int(metrics, \(x) if (is.null(x)) NA_integer_ else nrow(x)),
           metrics_has_depth = map_lgl(metrics, \(x) !is.null(x) && "depth" %in% colnames(x))) %>%
    select(-metrics, -file_names)
  rm(cm); invisible(gc(verbose = FALSE))

  flags <- flags %>%
    bind_cols(bind_rows(furrr::future_map(flags$afpd_dir, lf_read_rd)),
              bind_rows(furrr::future_map(flags$afpd_dir, lf_file_counts)) %>%
                select(raw_ranked)) %>%
    left_join(lf_bw_rows(), by = "p1_name")

  ledger_path <- fs::path(tree$work, "queue_ledger.tsv")
  ledger <- if (fs::file_exists(ledger_path)) {
    data.table::fread(ledger_path, sep = "\t") %>% as_tibble()
  } else tibble(afpd_dir_name = character(), attempts = integer(),
                last_job_dir = character(), last_submitted = character())
  flags <- flags %>%
    left_join(ledger %>% select(afpd_dir_name, attempts), by = "afpd_dir_name") %>%
    mutate(attempts = ifelse(is.na(attempts), 0L, as.integer(attempts)))

  ## ABSOLUTE completeness, cross-checked against ranking_debug.json, rather
  ## than the self-consistency of num_ark / num_con / num_pae. A short
  ## prediction set can never read as complete here.
  flags <- flags %>%
    mutate(af_ready = has_json_debug & rd_valid &
                      n_iptm == tree$n_pred & n_iptmptm == tree$n_pred &
                      n_order == tree$n_pred &
                      num_ark == tree$n_pred & num_con == tree$n_pred &
                      num_pae == tree$n_pred & num_xtr == tree$n_xtr &
                      num_models == tree$n_models,
           af_ready = af_ready %in% TRUE,
           capable  = !is.na(bw_rows) & bw_rows > 0,
           status = case_when(
             num_ark == 0 & raw_ranked > 0                    ~ "unrenamed_files",
             !capable                                          ~ "incapable_no_bw",
             !af_ready                                         ~ "incomplete_af",
             metrics_rows == tree$n_pred & metrics_has_depth    ~ "done",
             metrics_rows == tree$n_pred & !metrics_has_depth   ~ "incapable_no_depth",
             attempts >= max_attempts                          ~ "attempts_exhausted",
             TRUE                                              ~ "queue"))

  summary <- flags %>% count(status) %>% arrange(desc(n))
  print(summary)
  message("expected per directory: ", tree$n_pred,
          " predictions (ark = con = pae = length(iptm)), ", tree$n_xtr,
          " xtr, num_models ", tree$n_models)
  if (any(flags$status == "incomplete_af"))
    message("incomplete_af: short or unfinished prediction sets. Not queued and NOT marked\n",
            "done -- they return as 'queue' once the restart fills them in.")
  if (any(flags$status == "incapable_no_bw"))
    message("incapable_no_bw: receptor has no usable 'bw: full_table', so compute_RLdists()\n",
            "can never run. Reported and excluded, never re-queued.")
  if (any(flags$status == "incapable_no_depth"))
    message("incapable_no_depth: ATTEMPTED -- metrics_v2.csv exists with the full ",
            tree$n_pred, " rows\nbut no depth column (every row seq_match2 == 'different', or an empty bw\n",
            "table). Reported and excluded, never re-queued.")
  if (any(flags$status == "attempts_exhausted"))
    message("attempts_exhausted: queued ", max_attempts, "x with no usable metrics_v2.csv.\n",
            "Read the worker logs under the job dirs named in ", ledger_path)
  if (any(flags$status == "unrenamed_files"))
    message("unrenamed_files: run lf_resume_files(tree, lf_census(tree)) first.")

  lf_write_tsv(flags, fs::path(tree$work, paste0("gate_", lf_stamp(), ".tsv")))
  invisible(list(flags = flags, summary = summary, ledger_path = ledger_path,
                 max_attempts = max_attempts, tree_label = tree$label))
}

## === STEP 4  submit (sbatch is IRREVERSIBLE) ==============================

lf_submit <- function(tree, gate, group_size = NULL, dry_run = TRUE) {
  lf_hr("submit ", tree$label, if (dry_run) "  (dry run)" else "  IRREVERSIBLE")
  lf_assert(gate$tree_label == tree$label, "gate is from a different tree")
  ## extract_metrics.sh allows 12h and 16G per array task and the worker loops
  ## over its group serially, so budget by models rather than directories:
  ## ~50 models per task, i.e. 10 directories for A and C, 2 for B.
  if (is.null(group_size)) group_size <- max(1L, 50L %/% tree$n_pred)

  q <- gate$flags %>% filter(status == "queue") %>%
         select(afpd_dir, afpd_dir_name, p1_name, p2_name, run_dir)
  if (nrow(q) == 0) { message("nothing to queue"); return(invisible(NULL)) }

  lf_assert(fs::dir_exists(q$afpd_dir), "queued directory does not exist", q$afpd_dir)
  lf_assert(!is.na(q$p1_name) & nzchar(q$p1_name), "empty p1_name", q$afpd_dir_name)
  lf_assert(!is.na(q$p2_name) & nzchar(q$p2_name), "empty p2_name", q$afpd_dir_name)
  lf_assert(fs::path_file(q$run_dir) == tree$label,
            "run_dir would give submit_metric_jobs() the wrong run_id")
  lf_assert(!grepl("[\t\n]", paste0(q$afpd_dir, q$p1_name, q$p2_name, q$run_dir)),
            "a tab or newline would corrupt the TSV group file")
  lf_assert(anyDuplicated(q$afpd_dir) == 0, "duplicate directory in the queue")

  job_dir <- fs::path(LF_JOB_ROOT, paste0(tree$label, "_", lf_stamp()))
  lf_assert(!fs::dir_exists(job_dir), "job dir already exists", job_dir)

  message(nrow(q), " directories, group_size ", group_size, " -> ",
          ceiling(nrow(q) / group_size), " array tasks")
  submit_metric_jobs(runs = q, group_size = group_size, job_dir = job_dir, dry_run = dry_run)

  if (dry_run) {
    message("read ", job_dir, "/group_1.txt, then re-run with dry_run = FALSE")
    return(invisible(job_dir))
  }

  ## Record the attempt, so the next pass can tell "never attempted" from
  ## "attempted and failed".
  old <- if (fs::file_exists(gate$ledger_path)) {
    data.table::fread(gate$ledger_path, sep = "\t") %>% as_tibble()
  } else tibble(afpd_dir_name = character(), attempts = integer(),
                last_job_dir = character(), last_submitted = character())
  prev <- setNames(as.integer(old$attempts), old$afpd_dir_name)
  upd <- tibble(afpd_dir_name = q$afpd_dir_name,
                attempts = as.integer(ifelse(is.na(prev[q$afpd_dir_name]), 0L,
                                             prev[q$afpd_dir_name])) + 1L,
                last_job_dir = as.character(job_dir),
                last_submitted = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"))
  lf_write_tsv(bind_rows(old %>% filter(!afpd_dir_name %in% upd$afpd_dir_name), upd),
               gate$ledger_path)
  message("submitted; attempt ledger at ", gate$ledger_path)
  invisible(job_dir)
}
