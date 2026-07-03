# =============================================================================
# Methanotroph Taxonomy Reconciliation
# =============================================================================
#
# Maps equivalent taxonomy strings across databases (SILVA, NCBI/Kaiju, GTDB)
# to a single canonical name. Plug in new synonyms as you encounter them.
#
# Usage:
#   1. Edit the mapping table below to add/remove synonyms
#   2. Source this script
#   3. Call reconcile_taxonomy() on any phyloseq object or data frame
#
# =============================================================================

library(dplyr)
library(stringr)
library(tibble)

# =============================================================================
# MAPPING TABLE — edit this to add synonyms
# =============================================================================
# Each row: a database-specific string and the canonical name it should map to.
# canonical_name = what you want reported in your final figures/tables.
# db_name        = what actually appears in SILVA / Kaiju / GTDB output.
# rank           = which taxonomic rank this mapping applies to.
# db_source      = which database uses this string (for your reference;
#                  not used in matching — all db_names are checked against
#                  all inputs regardless of source).
#
# Add rows freely. If a name is already correct (canonical = db), include it
# so the table is a complete inventory of expected methanotroph names.

taxonomy_map <- tribble(
  ~rank,    ~db_source,   ~db_name,                          ~canonical_name,

  # =========================================================================
  # ORDER
  # =========================================================================
  "order",  "all",        "Methylococcales",                 "Methylococcales",

  # Rhizobiales/Hyphomicrobiales — contains Type II methanotroph families.
  # Rhizobiales was reclassified as Hyphomicrobiales (Oren & Garrity 2021)
  # but many databases still use Rhizobiales.
  "order",  "all",        "Rhizobiales",                    "Rhizobiales/Hyphomicrobiales",
  "order",  "NCBI",       "Hyphomicrobiales",               "Rhizobiales/Hyphomicrobiales",
  "order",  "SILVA",      "Hyphomicrobiales",               "Rhizobiales/Hyphomicrobiales",
  "order",  "GTDB",       "Hyphomicrobiales",               "Rhizobiales/Hyphomicrobiales",

  # Verrucomicrobial methanotrophs
  "order",  "SILVA",      "Methylacidiphilales",             "Methylacidiphilales",
  "order",  "NCBI",       "Methylacidiphilales",             "Methylacidiphilales",
  "order",  "GTDB",       "Methylacidiphilales",             "Methylacidiphilales",

  # Anaerobic methanotrophs (ANME)
  # ANME-2/3 are within Methanosarcinales; ANME-1 is Ca. Methanophagales
  # ANME-2/3 are within Methanosarcinales; ANME-1 is Ca. Methanophagales
  "order",  "NCBI",       "Methanosarcinales",               "Methanosarcinales",
  "order",  "SILVA",      "Methanosarcinales",               "Methanosarcinales",
  "order",  "SILVA",      "Methanosarciniales",              "Methanosarcinales",


  "order",  "NCBI",       "Methanophagales",                 "Ca. Methanophagales",
  "order",  "NCBI",       "Ca. Methanophagales",             "Ca. Methanophagales",
  "order",  "SILVA",      "ANME-1",                          "Ca. Methanophagales",
  "order",  "GTDB",       "Methanophagales",                 "Ca. Methanophagales",

  # NC10 — Candidatus Methylomirabilis
  # NCBI has sparse lineage: Bacteria;candidate division NC10;NA;NA;NA;Ca. Methylomirabilis
  # SILVA may use Methylomirabilota (phylum) or Methylomirabilia (class)
  "order",  "NCBI",       "Methylomirabilales",              "Methylomirabilales",
  "order",  "SILVA",      "Methylomirabilales",              "Methylomirabilales",
  "order",  "GTDB",       "Methylomirabilales",              "Methylomirabilales",

  # Phylum-level synonyms for NC10
  "phylum", "NCBI",       "candidate division NC10",         "Methylomirabilota",
  "phylum", "SILVA",      "Methylomirabilota",               "Methylomirabilota",
  "phylum", "GTDB",       "Methylomirabilota",               "Methylomirabilota",

  # =========================================================================
  # FAMILY — the problematic rank
  # =========================================================================
  # Methylococcaceae is the only validly published name (ICNP).
  # Methylomonadaceae (Leadbetter 1974) was never validly published and is
  # listed as a heterotypic synonym on LPSN. Some databases (SILVA, GTDB)
  # use it based on phylogenomic proposals (Orata et al. 2018).
  # Canonical name uses both separated by / to acknowledge the ambiguity.

  "family", "NCBI",       "Methylococcaceae",                "Methylococcaceae/Methylomonadaceae",
  "family", "SILVA",      "Methylococcaceae",                "Methylococcaceae/Methylomonadaceae",
  "family", "GTDB",       "Methylococcaceae",                "Methylococcaceae/Methylomonadaceae",
  "family", "NCBI",       "Methylomonadaceae",               "Methylococcaceae/Methylomonadaceae",
  "family", "SILVA",      "Methylomonadaceae",               "Methylococcaceae/Methylomonadaceae",
  "family", "GTDB",       "Methylomonadaceae",               "Methylococcaceae/Methylomonadaceae",

  "family", "all",        "Methylothermaceae",               "Methylothermaceae",
  "family", "all",        "Methylocystaceae",                "Methylocystaceae",
  "family", "all",        "Beijerinckiaceae",                "Beijerinckiaceae",
  "family", "all",        "Methylacidiphilaceae",            "Methylacidiphilaceae",
  "family", "all",        "Methanoperedenaceae",             "Methanoperedenaceae",
  "family", "NCBI",       "Candidatus Methanoperedenaceae",  "Methanoperedenaceae",

  # ANME families (Chadwick et al. 2022)
  "family", "all",        "Methanosarcinaceae",              "Methanosarcinaceae",
  "family", "all",        "Methanocomedenaceae",             "Methanocomedenaceae",
  "family", "all",        "Methanogasteraceae",              "Methanogasteraceae",
  "family", "all",        "Methanophagaceae",                "Methanophagaceae",
  "family", "all",        "Methylomirabilaceae",             "Methylomirabilaceae",

  # =========================================================================
  # GENUS — most stable rank; add synonyms / reclassifications here
  # =========================================================================

  # --- Methylococcaceae sensu stricto genera ---
  "genus",  "all",        "Methylococcus",                   "Methylococcus",
  "genus",  "all",        "Methylocaldum",                   "Methylocaldum",
  "genus",  "all",        "Methylothermus",                  "Methylothermus",
  "genus",  "all",        "Methylomagnum",                   "Methylomagnum",
  "genus",  "all",        "Methyloterrigena",                "Methyloterrigena",

  # --- Methylomonadaceae genera (may appear under Methylococcaceae in old dbs) ---
  "genus",  "all",        "Methylomonas",                    "Methylomonas",
  "genus",  "all",        "Methylobacter",                   "Methylobacter",
  "genus",  "all",        "Methylosarcina",                  "Methylosarcina",
  "genus",  "all",        "Methylovulum",                    "Methylovulum",
  "genus",  "all",        "Methylomarinum",                  "Methylomarinum",
  "genus",  "all",        "Methylomarinovum",                "Methylomarinovum",
  "genus",  "all",        "Methyloprofundus",                "Methyloprofundus",
  "genus",  "all",        "Methyloparacoccus",               "Methyloparacoccus",
  "genus",  "all",        "Methylogaea",                     "Methylogaea",
  "genus",  "all",        "Methylosoma",                     "Methylosoma",
  "genus",  "all",        "Methyloglobulus",                  "Methyloglobulus",
  "genus",  "all",        "Methylosphaera",                  "Methylosphaera",
  "genus",  "all",        "Methyloholbius",                  "Methyloholbius",

  # Methylomicrobium was split — some strains moved to Methylotuvimicrobium
  "genus",  "all",        "Methylomicrobium",                "Methylomicrobium",
  "genus",  "NCBI",       "Methylotuvimicrobium",            "Methylotuvimicrobium",
  "genus",  "GTDB",       "Methylotuvimicrobium",            "Methylotuvimicrobium",
  # If you want to lump them back together, change canonical to "Methylomicrobium"
  # "genus", "NCBI",      "Methylotuvimicrobium",            "Methylomicrobium",

  # Methylicorpusculum — split from Methylomonas
  "genus",  "NCBI",       "Methylicorpusculum",              "Methylicorpusculum",

  # --- Type II methanotrophs (Alphaproteobacteria) ---
  "genus",  "all",        "Methylosinus",                    "Methylosinus",
  "genus",  "all",        "Methylocystis",                   "Methylocystis",
  "genus",  "all",        "Methylocella",                    "Methylocella",
  "genus",  "all",        "Methylocapsa",                    "Methylocapsa",
  "genus",  "all",        "Methyloferula",                   "Methyloferula",

  # --- Verrucomicrobial methanotrophs ---
  "genus",  "all",        "Methylacidiphilum",               "Methylacidiphilum",
  "genus",  "all",        "Methylacidimicrobium",            "Methylacidimicrobium",

  # --- Anaerobic methanotrophs (ANME) ---
  # SILVA uses ANME cluster names; NCBI/Chadwick 2022 uses formal genus names
  # ANME-1 = Ca. Methanophagales (order) / Methanophagaceae (family)
  # ANME-2a/2b = Methanocomedenaceae (family), genus Methanocomedens / Methanomarinus
  # ANME-2c = Methanogasteraceae (family), genus Methanogaster
  # ANME-2d = Methanoperedenaceae (family), genus Methanoperedens
  # ANME-3 = Methanosarcinaceae (family), genus Methanovorans

  "genus",  "NCBI",       "Candidatus Methanoperedens",      "Ca. Methanoperedens",
  "genus",  "SILVA",      "Candidatus Methanoperedens",      "Ca. Methanoperedens",
  "genus",  "SILVA",      "ANME-2d",                         "Ca. Methanoperedens",
  "genus",  "GTDB",       "ANME-2d",                         "Ca. Methanoperedens",
  "genus",  "NCBI",       "Methanoperedens",                 "Ca. Methanoperedens",

  "genus",  "SILVA",      "ANME-2a-2b",                      "Methanocomedens/Methanomarinus",
  "genus",  "NCBI",       "Methanocomedens",                 "Methanocomedens",
  "genus",  "NCBI",       "Methanomarinus",                  "Methanomarinus",

  "genus",  "SILVA",      "ANME-2c",                         "Methanogaster",
  "genus",  "NCBI",       "Methanogaster",                   "Methanogaster",

  "genus",  "SILVA",      "ANME-3",                          "Methanovorans",
  "genus",  "NCBI",       "Methanovorans",                   "Methanovorans",

  # ANME-1 — SILVA splits into ANME-1a, ANME-1b at class level with
  # Incertae Sedis at family/genus. Formal names: order Methanophagales,
  # family Methanophagaceae, genera Methanophaga / Methanoalium
  "genus",  "SILVA",      "ANME-1",                          "Methanophaga",
  "genus",  "SILVA",      "ANME-1a",                         "Methanophaga",
  "genus",  "SILVA",      "ANME-1b",                         "Methanophaga",
  "genus",  "NCBI",       "Methanophaga",                    "Methanophaga",
  "genus",  "NCBI",       "Methanoalium",                    "Methanoalium",

  # --- NC10 ---
  "genus",  "NCBI",       "Candidatus Methylomirabilis",     "Ca. Methylomirabilis",
  "genus",  "SILVA",      "Candidatus Methylomirabilis",     "Ca. Methylomirabilis",
  "genus",  "SILVA",      "Methylomirabilis",                "Ca. Methylomirabilis",
  "genus",  "GTDB",       "Methylomirabilis",                "Ca. Methylomirabilis"
)


# =============================================================================
# GENUS-TO-FAMILY CORRECTION TABLE
# =============================================================================
# Use this to fix family assignments when the database has stale lineages
# (e.g., Kaiju putting Methylomonas under Methylococcaceae).
# This is applied AFTER the genus name is reconciled.

genus_to_correct_lineage <- tribble(
  ~canonical_genus,         ~correct_order,                     ~correct_family,
  # --- Type I methanotrophs (Gammaproteobacteria) ---
  "Methylomonas",           "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylobacter",          "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylomicrobium",       "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylotuvimicrobium",   "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylosarcina",         "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylovulum",           "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylomarinum",         "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylomarinovum",       "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methyloprofundus",       "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methyloparacoccus",      "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylogaea",            "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylosoma",            "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methyloglobulus",        "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylosphaera",         "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methyloholbius",         "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylicorpusculum",     "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylococcus",          "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylocaldum",          "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methylothermus",         "Methylococcales",                  "Methylothermaceae",
  "Methylomagnum",          "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",
  "Methyloterrigena",       "Methylococcales",                  "Methylococcaceae/Methylomonadaceae",

  # --- Methylothermaceae genera (separate family within Methylococcales) ---
  "Methylohalobius",        "Methylococcales",                  "Methylothermaceae",
  "Methylomarinovum",       "Methylococcales",                  "Methylothermaceae",

  # --- Type II methanotrophs (Alphaproteobacteria) ---
  "Methylosinus",           "Rhizobiales/Hyphomicrobiales",     "Methylocystaceae",
  "Methylocystis",          "Rhizobiales/Hyphomicrobiales",     "Methylocystaceae",
  "Methylocella",           "Rhizobiales/Hyphomicrobiales",     "Beijerinckiaceae",
  "Methylocapsa",           "Rhizobiales/Hyphomicrobiales",     "Beijerinckiaceae",
  "Methyloferula",          "Rhizobiales/Hyphomicrobiales",     "Beijerinckiaceae",

  # --- Verrucomicrobial methanotrophs ---
  "Methylacidiphilum",      "Methylacidiphilales",              "Methylacidiphilaceae",
  "Methylacidimicrobium",   "Methylacidiphilales",              "Methylacidiphilaceae",

  # --- NC10 (Candidatus Methylomirabilis) ---
  # Kaiju/NCBI lineage has NA at class/order/family; fill them in
  "Ca. Methylomirabilis",   "Methylomirabilales",               "Methylomirabilaceae",

  # --- ANME-1 (Ca. Methanophagales) ---
  "Methanophaga",           "Ca. Methanophagales",              "Methanophagaceae",
  "Methanoalium",           "Ca. Methanophagales",              "Methanophagaceae",

  # --- ANME-2a (Methanosarcinales) ---
  "Methanocomedens",        "Methanosarcinales",                "Methanocomedenaceae",

  # --- ANME-2b (Methanosarcinales) ---
  "Methanomarinus",         "Methanosarcinales",                "Methanocomedenaceae",

  # --- SILVA lumps 2a+2b together; this catches the slash canonical name ---
  "Methanocomedens/Methanomarinus", "Methanosarcinales",        "Methanocomedenaceae",

  # --- ANME-2c (Methanosarcinales) ---
  "Methanogaster",          "Methanosarcinales",                "Methanogasteraceae",

  # --- ANME-2d (Methanosarcinales) ---
  "Ca. Methanoperedens",    "Methanosarcinales",                "Methanoperedenaceae",

  # --- ANME-3 (Methanosarcinales) ---
  "Methanovorans",          "Methanosarcinales",                "Methanosarcinaceae"
)


# =============================================================================
# RECONCILIATION FUNCTIONS
# =============================================================================

#' Reconcile a single taxonomy column in a data frame
#'
#' @param df        Data frame with a taxonomy column
#' @param col       Name of the column to reconcile (quoted string)
#' @param rank      Rank of that column: "order", "family", or "genus"
#' @param map       The mapping table (defaults to taxonomy_map)
#' @return          The data frame with the column values replaced by canonical names
reconcile_column <- function(df, col, rank,
                             map = taxonomy_map) {
  rank_map <- map %>%
    filter(rank == !!rank) %>%
    select(db_name, canonical_name) %>%
    distinct()

  # Build a named lookup vector: db_name -> canonical_name
  lookup <- setNames(rank_map$canonical_name, rank_map$db_name)

  df[[col]] <- ifelse(
    df[[col]] %in% names(lookup),
    lookup[df[[col]]],
    df[[col]]  # keep as-is if not in map
  )
  df
}


#' Reconcile taxonomy in a phyloseq tax_table
#'
#' @param physeq    A phyloseq object
#' @param ranks     Named list mapping rank name -> column name in tax_table.
#'                  Defaults assume standard column names.
#' @param map       The mapping table
#' @return          The phyloseq object with reconciled taxonomy
reconcile_phyloseq <- function(physeq,
                               ranks = NULL,
                               map = taxonomy_map) {
  tax_df <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)

  # Auto-detect rank columns
  if (is.null(ranks)) {
    all_cols <- colnames(tax_df)
    ranks <- list()
    if ("order"  %in% all_cols) ranks[["order"]]  <- "order"
    if ("family" %in% all_cols) ranks[["family"]] <- "family"
    if ("genus"  %in% all_cols) ranks[["genus"]]  <- "genus"
  }

  for (rank_name in names(ranks)) {
    col_name <- ranks[[rank_name]]
    if (col_name %in% colnames(tax_df)) {
      tax_df <- reconcile_column(tax_df, col_name, rank_name, map)
    }
  }

  tax_table(physeq) <- tax_table(as.matrix(tax_df))
  physeq
}


#' Fix order and family assignments based on genus
#'
#' Fills in missing or incorrect order/family using the genus_to_correct_lineage table.
#' Especially important for lineages like NC10 where Kaiju/NCBI leaves
#' intermediate ranks as NA/Unclassified.
#'
#' @param physeq         A phyloseq object (genus column should already be reconciled)
#' @param genus_col      Name of the genus column
#' @param order_col      Name of the order column
#' @param family_col     Name of the family column
#' @param corrections    The correction table
#' @param overwrite      If TRUE, overwrite all order/family values from the table.
#'                       If FALSE (default), only fill in "Unclassified" or NA values.
#' @return               Updated phyloseq
fix_lineage_from_genus <- function(physeq,
                                   genus_col = "genus",
                                   order_col = "order",
                                   family_col = "family",
                                   corrections = genus_to_correct_lineage,
                                   overwrite = FALSE) {
  tax_df <- as.data.frame(tax_table(physeq), stringsAsFactors = FALSE)

  missing_cols <- setdiff(c(genus_col, order_col, family_col), colnames(tax_df))
  if (length(missing_cols) > 0) {
    warning("Column(s) not found: ", paste(missing_cols, collapse = ", "),
            ". Skipping lineage correction.")
    return(physeq)
  }

  order_lookup  <- setNames(corrections$correct_order,  corrections$canonical_genus)
  family_lookup <- setNames(corrections$correct_family, corrections$canonical_genus)

  has_genus <- tax_df[[genus_col]] %in% names(order_lookup)

  if (overwrite) {
    # Overwrite all matching rows
    tax_df[[order_col]][has_genus]  <- order_lookup[tax_df[[genus_col]][has_genus]]
    tax_df[[family_col]][has_genus] <- family_lookup[tax_df[[genus_col]][has_genus]]
  } else {
    # Only fill where order/family is Unclassified or NA
    needs_order  <- has_genus & (tax_df[[order_col]]  %in% c("Unclassified", NA, "NA"))
    needs_family <- has_genus & (tax_df[[family_col]] %in% c("Unclassified", NA, "NA"))

    tax_df[[order_col]][needs_order]   <- order_lookup[tax_df[[genus_col]][needs_order]]
    tax_df[[family_col]][needs_family] <- family_lookup[tax_df[[genus_col]][needs_family]]
  }

  tax_table(physeq) <- tax_table(as.matrix(tax_df))
  physeq
}


#' Quick summary: show what got remapped
#'
#' @param physeq_before   phyloseq before reconciliation
#' @param physeq_after    phyloseq after reconciliation
#' @param rank            which rank to compare
compare_reconciliation <- function(physeq_before, physeq_after, rank) {
  before <- as.data.frame(tax_table(physeq_before))[[rank]]
  after  <- as.data.frame(tax_table(physeq_after))[[rank]]

  changed <- before != after & !is.na(before) & !is.na(after)

  if (sum(changed) == 0) {
    message("No changes at ", rank, " level.")
    return(invisible(NULL))
  }

  changes <- tibble(
    original    = before[changed],
    reconciled  = after[changed]
  ) %>%
    distinct() %>%
    arrange(original)

  message(sum(changed), " entries remapped at ", rank, " level:")
  print(changes)
  invisible(changes)
}


# =============================================================================
# EXAMPLE USAGE
# =============================================================================
#
# # After building your phyloseq objects (physeq_16s, physeq_mg):
#
# # 1. Reconcile genus/order/family names to canonical forms
# physeq_16s <- reconcile_phyloseq(physeq_16s)
# physeq_mg  <- reconcile_phyloseq(physeq_mg)
#
# # 2. Fill in missing order/family from genus (especially for Kaiju/NC10)
# #    By default only fills "Unclassified" or NA values.
# #    Set overwrite = TRUE to force all order/family to match the table.
# physeq_mg  <- fix_lineage_from_genus(physeq_mg)
# physeq_16s <- fix_lineage_from_genus(physeq_16s)
#
# # 3. Check what changed
# compare_reconciliation(physeq_mg_original, physeq_mg, "order")
# compare_reconciliation(physeq_mg_original, physeq_mg, "family")
# compare_reconciliation(physeq_mg_original, physeq_mg, "genus")
#
# # 4. To add a new synonym you encounter:
# #    Just add a row to taxonomy_map above, e.g.:
# #    "genus", "GTDB", "Methylomonas_A", "Methylomonas",
# #    Then re-source this script and re-run reconcile_phyloseq().
#
# # 5. To add a genus with missing higher ranks:
# #    Add a row to genus_to_correct_lineage, e.g.:
# #    "NewGenus", "SomeOrder", "SomeFamily",
# #    Then re-run fix_lineage_from_genus().
