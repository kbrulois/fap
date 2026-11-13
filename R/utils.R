

residue_pairs <- list(
  favorable = list(
    hydrogen_bonds = list(
      c("D", "R"), c("E", "K"), c("S", "D"), c("Y", "H"),
      c("S", "E"), c("T", "E"), c("S", "R"), c("N", "Q"),
      c("H", "D"), c("H", "E")
    ),
    salt_bridges = list(
      c("D", "K"), c("E", "R"), c("D", "H"), c("E", "K"),
      c("D", "R"), c("E", "H")
    ),
    hydrophobic_interactions = combn(
      c("A", "V", "L", "I", "M", "F", "W", "Y"), 2, simplify = FALSE
    ),
    aromatic_stacking = list(
      c("F", "Y"), c("Y", "W"), c("F", "W")
    ),
    cation_pi = list(
      c("F", "R"), c("F", "K"), c("Y", "R"), c("Y", "K"), c("W", "R"), c("W", "K")
    )
  ),
  unfavorable = list(
    steric_clashes = list(
      c("W", "W"), c("F", "W"), c("F", "F"), c("I", "L"), c("L", "L"), c("Y", "W")
    ),
    electrostatic_repulsion = list(
      c("D", "E"), c("K", "R"), c("D", "D"), c("E", "E"), c("K", "K"), c("R", "R")
    ),
    desolvation_penalty = list(
      c("S", "S"), c("T", "T"), c("N", "N"), c("Q", "Q")
    )
  )
)

residue_pairs <- c(residue_pairs$favorable, residue_pairs$unfavorable)

res_pairs <- lapply(names(residue_pairs), \(x) tibble(pairs = map(residue_pairs[[x]], \(x) c(paste0(x[1], x[2]), paste0(x[2], x[1]))) %>% do.call(args = ., what = base::c),
                                         name = x)) %>%
  bind_rows


alpha_fold_AA_order <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'X', '-')

#from fig1
cp_bw <- c("1.28", "1.35", "1.39", "2.53", "2.56", "2.57", "2.59", "2.60", "2.61", "2.63", "2.64", "3.21", "3.26", "3.28", "3.29", "3.30", "3.32", "3.33", "3.34", "3.36", "3.37", "3.40", "4.56", "4.57", "4.58", "4.59", "4.60", "4.61", "5.33", "5.35", "5.36", "5.37", "5.39", "5.40", "5.43", "5.44", "5.46", "5.461", "5.47", "6.44", "6.48", "6.51", "6.52", "6.54", "6.55", "6.58", "6.59", "6.61", "6.62", "7.24", "7.27", "7.30", "7.31", "7.33", "7.34", "7.35", "7.37", "7.38", "7.39", "7.41", "7.42")

residue_pairs <- list(

  # Favorable interactions: Hydrogen bonds, salt bridges, hydrophobic packing
  favorable = list(
    hydrogen_bonds = list(
      c("D", "R"), c("E", "K"), c("S", "D"), c("Y", "H")
    ),
    salt_bridges = list(
      c("D", "K"), c("E", "R"), c("D", "H"), c("E", "K")
    ),
    hydrophobic_interactions = list(
      c("L", "I"), c("F", "L"), c("V", "M"), c("W", "A")
    ),
    aromatic_stacking = list(
      c("F", "Y"), c("Y", "W"), c("F", "W")
    )
  ),

  # Unfavorable interactions: Steric clashes, repulsive electrostatics, desolvation penalties
  unfavorable = list(
    steric_clashes = list(
      c("W", "W"), c("F", "W"), c("I", "I"), c("L", "L")
    ),
    electrostatic_repulsion = list(
      c("D", "E"), c("K", "R"), c("D", "D"), c("E", "E")
    ),
    desolvation_penalty = list(
      c("S", "S"), c("T", "T"), c("D", "D")
    )
  )
)

residue_pairs <- lapply(residue_pairs, \(x) unname(do.call(c, x)))


residue_pairs <- list(

  # Favorable interactions: Hydrogen bonds, salt bridges, hydrophobic packing
  favorable = list(
    hydrogen_bonds = list(
      c("D", "R"), c("E", "K"), c("S", "D"), c("Y", "H")
    ),
    salt_bridges = list(
      c("D", "K"), c("E", "R"), c("D", "H"), c("E", "K")
    ),
    hydrophobic_interactions = list(
      c("L", "I"), c("F", "L"), c("V", "M"), c("W", "A")
    ),
    aromatic_stacking = list(
      c("F", "Y"), c("Y", "W"), c("F", "W")
    )
  ),

  # Unfavorable interactions: Steric clashes, repulsive electrostatics, desolvation penalties
  unfavorable = list(
    steric_clashes = list(
      c("W", "W"), c("F", "W"), c("I", "I"), c("L", "L")
    ),
    electrostatic_repulsion = list(
      c("D", "E"), c("K", "R"), c("D", "D"), c("E", "E")
    ),
    desolvation_penalty = list(
      c("S", "S"), c("T", "T"), c("D", "D")
    )
  )
)


bw_segs <- tibble::tibble(protein_segment = c("N-term", "TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4",
  "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "C-term"))

dist_colors <- c(setNames("#D2AF81FF", "0-2A"),
                 setNames("#FD7446FF", "2-4A"),
                 setNames("#197EC0FF", "4-6A"),
                 setNames("#FED439FF", "6-8A"),
                 setNames("#8A919780", "> 8A"))

cons_colors <- c(setNames("#C80813FF", "< 1"),
                 setNames("#FD7446FF", "1-2"),
                 setNames("#FD8CC1FF", "2-3"),
                 setNames("#46732E80", "3-4"),
                 setNames("#197EC0FF", "> 4"))

unique_non_na <- function(x) unique(x[!is.na(x)])

mod_codes <- function() {
  expand.grid(letters, letters) %>%
                dplyr::mutate(codes = paste0(Var2, Var1)) %>%
                dplyr::pull(codes)
}

mod_codes <- mod_codes()

rlx_conv <- c(setNames("u", "unrelaxed"),
              setNames("r", "relaxed"))

file_type_conv <- c(setNames("pae", "pae"),
                    setNames("con", "confidence"),
                    setNames("ark", "ranked"),
                    setNames("xtr", "unrelaxed"),
                    setNames("xtr", "relaxed"),
                    setNames("res", "result"),
                    setNames("plt", "pae_plot_ranked"),
                    setNames("spc", "paeSPOC"))

file_type_conv_af3 <- c(setNames("pae", "confidences"),
                        setNames("con", "summary_confidences"),
                        setNames("ark", "model"))

aa_df <- data.frame(
  OneLetter = c(
    "A", "R", "N", "D", "C",
    "E", "Q", "G", "H", "I",
    "L", "K", "M", "F", "P",
    "S", "T", "W", "Y", "V",
    "U", "O", "X", "B", "Z", "J"
  ),
  ThreeLetter = c(
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "XAA", "ASX", "GLX", "XLE"
  ),
  stringsAsFactors = FALSE
)


