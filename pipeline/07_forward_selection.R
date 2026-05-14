#!/usr/bin/env Rscript
# ============================================================================
# Forward selection on conditional MAGMA cell-type results.
#
# Algorithm from Duncan et al. 2025 (Nat Neurosci):
#   duncan_repo/MAGMA/5.forward_selection_condition_results.Rmd
# (R code authored by Tayden Li, Stanford)
#
# Modifications from Duncan's original:
#   1. File paths come from command-line args.
#   2. read_table uses comment="#" instead of skip=4 to robustly handle MAGMA
#      .gsa.out files with 3 or 4 # header lines (depending on whether the
#      GWAS has variable per-SNP N).
#   3. Adds a defensive swap of a/b columns to enforce the P_marg_a <= P_marg_b
#      invariant that Duncan's algorithm assumes. MAGMA's joint-pairs output
#      orders pair-rows by spec-matrix column position (lexicographic), NOT
#      by marginal significance, so naive alternating a/b assignment doesn't
#      always satisfy the invariant. Swap fixes this; algorithm semantics
#      unchanged ('a' = more marginally significant of the pair).
#   4. Independent cluster list written to file in addition to stdout.
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop(paste0("Usage: Rscript 07_forward_selection.R ",
                "<marginal.gsa.out> <joint_pairs.gsa.out> <output_indep.txt>"))
}
filename_marg <- args[1]
filename_cond <- args[2]
output_file   <- args[3]

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
})

# ---- Read inputs ----
df_cond <- read_table(filename_cond, comment="#", show_col_types=FALSE) %>%
    subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>%
    rename(P_cond = P)
df_marg <- read_table(filename_marg, comment="#", show_col_types=FALSE) %>%
    subset(select=-c(TYPE, NGENES, BETA, BETA_STD, SE)) %>%
    rename(P_marg = P) %>%
    arrange(P_marg)
n_sig <- ceiling(sqrt(nrow(df_cond)))
n_clusters <- nrow(df_marg)

# ---- Combine marginal and conditional ----
df_comb <- left_join(df_cond, df_marg, by = "VARIABLE")
df_comb <- df_comb %>%
    mutate(VarCode = rep(c("a", "b"), times=nrow(df_comb)/2)) %>%
    pivot_wider(names_from = VarCode,
                values_from = c(VARIABLE, P_cond, P_marg))

# ---- ADAPTATION: enforce P_marg_a <= P_marg_b invariant ----
needs_swap <- df_comb$P_marg_a > df_comb$P_marg_b
n_swapped <- sum(needs_swap)
if (n_swapped > 0) {
    cat(sprintf("Note: swapping a/b in %d / %d pairs to enforce P_marg_a <= P_marg_b\n",
                n_swapped, nrow(df_comb)))
    for (col_root in c("VARIABLE", "P_cond", "P_marg")) {
        a_col <- paste0(col_root, "_a")
        b_col <- paste0(col_root, "_b")
        tmp <- df_comb[[a_col]][needs_swap]
        df_comb[[a_col]][needs_swap] <- df_comb[[b_col]][needs_swap]
        df_comb[[b_col]][needs_swap] <- tmp
    }
}

# ---- Compute PS scores (after swap) ----
df_comb <- df_comb %>%
    mutate(PS_a = log10(P_cond_a) / log10(P_marg_a),
           PS_b = log10(P_cond_b) / log10(P_marg_b))
stopifnot(df_comb$P_marg_a <= df_comb$P_marg_b)

# ---- Pre-discard "dominated" clusters (verbatim from Duncan) ----
select_list  <- df_marg$VARIABLE[1:n_sig]
reverse_list <- df_comb$VARIABLE_a[df_comb$PS_a < 0.2 & df_comb$PS_b >= 0.2]
select_list  <- select_list[select_list %in% reverse_list == FALSE][-1]

# ---- Forward selection (verbatim from Duncan's Rmd) ----
indep_list <- c(toString(df_marg[1, 1]))
for (cur_var in select_list) {
    for (indep_var in indep_list) {
        df_cur <- df_comb[df_comb$VARIABLE_a == indep_var &
                          df_comb$VARIABLE_b == cur_var, ]
        if ((df_cur$PS_a >= 0.8 & df_cur$PS_b >= 0.8) |
            (df_cur$PS_a >= 0.5 & df_cur$PS_b >= 0.5 & df_cur$P_cond_b < 0.05)) {
            if (tail(indep_list, 1) == indep_var) {
                indep_list <- c(indep_list, cur_var)
            }
        } else break
    }
}

# ---- Write outputs ----
cat(indep_list, sep="\n", file=output_file)
cat(sprintf("\nIndependent clusters (%d):\n", length(indep_list)))
cat(paste(indep_list, collapse=", "), "\n")
