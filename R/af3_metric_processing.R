

clean_af3_file_names <- function(run_dir) {

  tmp <- tibble(files = fs::dir_ls(run_dir) %>% basename(),
                file_parts = map(files, ~stringr::str_split(., "_", simplify = TRUE))) %>%
    mutate(file_part_len = map_int(file_parts, length)) %>%
    mutate(parse_proteins(files, delim_proteins = "_", delim_ranges = "x", delim_start_end = "x")) %>%
    mutate(data_files = furrr::future_map(files, ~fs::dir_ls(.) %>% basename()))


  tmp <- tmp %>%
          mutate(num_files = map_int(data_files, length))


    mutate(new = case_when(length(file_parts) == 2, \(x) {
      if(length(x) == 2) return(tibble(file_name_new = stringr::str_c(x[,1:2], collapse = "_"),
                                       date = "0"))
      if(length(x) == 4) return(tibble(file_name_new = stringr::str_c(x[,1:2], collapse = "_"),
                                       date = stringr::str_c(x[,3:4], collapse = "")))
    })) %>%

    mutate(new = map(file_parts, \(x) {
      if(length(x) == 2) return(tibble(file_name_new = stringr::str_c(x[,1:2], collapse = "_"),
                                       date = "0"))
      if(length(x) == 4) return(tibble(file_name_new = stringr::str_c(x[,1:2], collapse = "_"),
                                       date = stringr::str_c(x[,3:4], collapse = "")))
    })) %>%
    unnest(new) %>%
    group_by(file_name_new) %>%
    mutate(keep = case_when(date == max(date) ~ "keep",
                            TRUE ~ "remove")) %>%
    ungroup %>%
    split(., .[["keep"]])


  if("remove" %in% names(tmp)) {
  tmp[["remove"]]  %>%
    mutate(unlink(paste0(input_path_models, "/", files), recursive = TRUE))
  }

  tmp[["keep"]] %>%
    filter(date != "0") %>%
    rename_dir(run_dir = run_dir,
               input = .,
               from = "files",
               to = "file_name_new")

}


import_raw_metrics_af3 <- function(dir_name,
                                   run_name = run_id,
                                   algorithm = "AF3",
                                   pdb = "_ark_.*.pdb$",
                                   pae = "_pae_.*.json$",
                                   code_model_rank = "_[A-Z]{4}[a-z]{7}_s\\d+m\\d+p\\d+_r\\d+") {

  message("importing raw metrics for ", dir_name)

  files <- list.files(paste(input_path_models, dir_name, sep = "/"))

  models <- tibble(pdb_files = grep(pdb, files, value = TRUE),
                   pae_files = grep(pae, files, value = TRUE),
                   code_model_rank = stringr::str_extract(pdb_files, code_model_rank),
                   rlx = stringr::str_extract(pdb_files, "_[ru]_") %>%
                     stringr::str_remove(., "^_") %>%
                     stringr::str_remove(., "_$"),
                   model_c = stringr::str_extract(code_model_rank, "s\\d+m\\d+p\\d+"),
                   rank = stringr::str_extract(code_model_rank, "_r\\d+$") %>% stringr::str_remove(., "^_r"),
                   model_num = stringr::str_extract(model_c, "m\\d+") %>% stringr::str_remove(., "^m"),
                   pred_num = stringr::str_extract(model_c, "p\\d+") %>% stringr::str_remove(., "^p"),
                   model_e = paste0("model_", model_num, "_multimer_v3_", "pred_", pred_num)) %>%
    mutate(code = stringr::str_extract(code_model_rank, "[A-Z]{4}[a-z]{7}"), .before = "pdb_files") %>%
    mutate(algorithm = algorithm, .before = "pdb_files") %>%
    mutate(run_name = run_name, .before = "pdb_files") %>%
    select(-code_model_rank)


  dat <- jsonlite::read_json(paste(input_path_models, dir_name, "ranking_debug.json", sep = "/"))

  dat <- tibble(model_e = names(dat[["iptm"]]),
                iptm = unlist(dat[["iptm"]]),
                `iptm+ptm` = unlist(dat[["iptm+ptm"]]))

  if(sum(models[["model_e"]] %in% dat[["model_e"]]) != nrow(models)) warning("file name models do not match ranking_debug.json")

  pairing_info <- parse_dirname(pairing_dir = dir_name,
                                delim_proteins = "_",
                                delim_ranges = "x",
                                delim_start_end = "x") %>%
    mutate(parsed_pair = map(parsed_pair,
                             ~pivot_wider(.,
                                          names_from=c("protein", "annotation"),
                                          values_from=value))) %>%
    unnest(parsed_pair)

  dat <- bind_cols(dat, pairing_info)

  dat <- left_join(models, dat, by = "model_e")

  dat %>%
    mutate(pdb_files = setNames(models[["pdb_files"]], models[["model_e"]])[model_e]) %>%

    mutate(pdb = map(pdb_files, ~bio3d::read.pdb(paste(input_path_models, dir_name, ., sep = "/")))) %>%

    mutate(pdb.xyz = map(pdb, parse_pdb)) %>%

    mutate(pae = map(setNames(models[["pae_files"]], models[["model_e"]])[model_e],
                     \(x) {
                       tmp <- jsonlite::read_json(paste(input_path_models, dir_name, x, sep = "/"))
                       tmp[[1]][["predicted_aligned_error"]]
                     }))

}


