# Shared constants for the cycdb functional gene scripts (9.1-9.3).
#
# Requires `project_root` to already be defined in the calling environment
# (each driver script sets it from its own location before sourcing this file).

# sample_id_map (raw McM2X_..._S## -> display sample name) is the single
# source of truth defined in step_6; reuse it here rather than redefining it.
source(file.path(project_root, "scripts", "step_6_metagenomic_taxonomy", "constants.R"))

diamond_dir         <- file.path(project_root, "cycdb_diamond_output")
gene_map_path        <- file.path(diamond_dir, "combined_id2genemap.csv")
exp_metadata_path    <- file.path(diamond_dir, "ant_exp_map_MA.csv")
pathway_counts_path  <- file.path(diamond_dir, "pathway_counts.rds")

# DESeq2 differential-testing figures (lollipop plots) go here rather than
# the shared figures/ folder, which holds VST-normalized visualizations
# (heatmaps) and the step_6/6 kaiju figures.
deseq_output_dir <- file.path(project_root, "deseq_output")
dir.create(deseq_output_dir, showWarnings = FALSE)

e_value_threshold  <- 1e-5
identity_threshold <- 60

diamond_column_names <- c(
  "query", "protein", "identity", "length", "mismatches", "gap_opens",
  "q_start", "q_end", "s_start", "s_end", "e_value", "bit_score",
  "qlen", "slen", "stitle", "pident_full"
)

# Gene -> pathway lookup, grouped by cycle.
lop <- list(

  # methane cycles
  cenmetpat = c('mtrA', 'mtrB', 'mtrC', 'mtrD', 'mtrE', 'mtrF', 'mtrG', 'mtrH',
                'mcrA', 'mcrB', 'mcrC', 'mcrD', 'mcrG', 'hdrA1', 'hdrB1', 'hdrC1',
                'hdrA2', 'hdrB2', 'hdrC2', 'hdrD', 'hdrE', 'mvhA', 'mvhD', 'mvhG',
                'ehaA', 'ehaB', 'ehaC', 'ehaD', 'ehaE', 'ehaF', 'ehaG', 'ehaH',
                'ehaI', 'ehaJ', 'ehaK', 'ehaL', 'ehaM', 'ehaN', 'ehaO', 'ehaP',
                'ehaQ', 'ehaR', 'ehbA', 'ehbB', 'ehbC', 'ehbD', 'ehbE', 'ehbF',
                'ehbG', 'ehbH', 'ehbI', 'ehbJ', 'ehbK', 'ehbL', 'ehbM', 'ehbN',
                'ehbO', 'ehbP', 'ehbQ', 'mbhL', 'mbhK', 'mbhJ', 'rnfA', 'rnfB',
                'rnfC', 'rnfD', 'rnfE', 'rnfG', 'echA', 'echB', 'echC', 'echD',
                'echE', 'echF', 'vhoA', 'vhoC', 'vhoG', 'vhtA', 'vhtC', 'vhtG',
                'vhuA', 'vhuU', 'vhuD', 'vhuG', 'vhcA', 'vhcD', 'vhcG', 'fpoA',
                'fpoB', 'fpoC', 'fpoD', 'fpoF', 'fpoH', 'fpoI', 'fpoJ', 'fpoK',
                'fpoL', 'fpoM', 'fpoN', 'fpoO', 'fqoA', 'fqoD', 'fqoF', 'fqoH',
                'fqoJ', 'fqoK', 'fqoL', 'fqoM', 'fqoN', 'frhA', 'frhB', 'frhD',
                'frhG', 'fruA', 'fruB', 'fruD', 'fruG', 'frcA', 'frcB', 'frcD',
                'frcG'),
  hydro_met = c('fmdA', 'fmdB', 'fmdC', 'fmdD', 'fmdE', 'fmdF', 'fwdA', 'fwdB',
                'fwdC', 'fwdD', 'fwdE', 'fwdF', 'fwdG', 'fwdH', 'ftr', 'mch',
                'mtdA', 'mtdB', 'hmd', 'mer', 'metF'),
  aceti_met = c('acs', 'acsA', 'acsC', 'acsD', 'ackA', 'pta', 'cdhC', 'cdhD',
                'cdhE', 'cdhA', 'cdhB', 'acdA', 'acdB', 'cooS', 'cooF'),
  methyl_met = c('mtaA', 'mtaB', 'mtaC', 'mtbA', 'mtbB', 'mtbC', 'mtsA', 'mtsB',
                 'torA', 'torC', 'torD', 'torY', 'torZ', 'mttB', 'mttC', 'mtmB', 'mtmC'),
  aom = c('fmdA', 'fmdB', 'fmdC', 'fmdD', 'fmdE', 'fmdF', 'fwdA', 'fwdB', 'fwdC',
          'fwdD', 'fwdE', 'fwdF', 'fwdG', 'fwdH', 'ftr', 'mch', 'mtdA', 'mtdB',
          'hmd', 'mer', 'metF', 'mtrA', 'mtrB', 'mtrC', 'mtrD', 'mtrE', 'mtrF',
          'mtrG', 'mtrH', 'mcrA', 'mcrB', 'mcrC', 'mcrD', 'mcrG', 'hdrA1', 'hdrB1',
          'hdrC1', 'hdrA2', 'hdrB2', 'hdrC2', 'hdrD', 'hdrE', 'rnfA', 'rnfB',
          'rnfC', 'rnfD', 'rnfE', 'rnfG', 'frhA', 'frhB', 'frhD', 'frhG', 'fruA',
          'fruB', 'fruD', 'fruG', 'frcA', 'frcB', 'frcD', 'frcG', 'fpoA',
          'fpoB', 'fpoC', 'fpoD', 'fpoF', 'fpoH', 'fpoI', 'fpoJ', 'fpoK',
          'fpoL', 'fpoM', 'fpoN', 'fpoO', 'fqoA', 'fqoD', 'fqoF', 'fqoH',
          'fqoJ', 'fqoK', 'fqoL', 'fqoM', 'fqoN', 'narG', 'narZ', 'narH',
          'narY', 'narB', 'napH', 'nrfH', 'nrfA', 'nxrA', 'nod', 'nirK',
          'nirS', 'cytC'),
  oxid_met_c1 = c('pqqA', 'pqqB', 'pqqC', 'pqqD', 'pqqE', 'pqqF',
                  'mmoX', 'mmoY', 'mmoZ', 'mmoB', 'mmoC', 'mmoD',
                  'pmoA', 'pmoB', 'pmoC', 'amoA', 'amoB', 'amoC',
                  'xoxF1', 'xoxF2', 'xoxF4', 'xoxF5', 'mxaF', 'mxaI', 'mxaJ',
                  'mxaG', 'mxaA', 'mxaC', 'mxaK', 'mxaL', 'mxaD', 'mauA', 'mauB',
                  'mauC', 'mauD', 'mauE', 'mauF', 'mgsA', 'mgsB', 'mgsC', 'mgdA',
                  'mgdB', 'mgdD', 'qhpA', 'dcmR', 'dcmA', 'mdh', 'tmm'),
  oxid_formaldehyde = c('gfa', 'frmA', 'frmB', 'fghA', 'fae', 'mtdB', 'mch',
                        'mdo', 'fdm', 'fdhA-K00148'),
  oxid_formate = c('fdwA', 'fdwB', 'fdwE', 'fdsB', 'fdsD', 'fdsG', 'fdoG',
                   'fdoH', 'fdoI', 'fdhF', 'fdhA', 'fdhB'),
  serine = c('fchA', 'mtdA', 'dfrA1', 'dfrA12', 'dfrA10', 'dfrA19', 'folA',
             'glyA', 'hprA', 'gck', 'eno', 'ppc', 'mtkA', 'mtkB', 'mcl', 'gpmI',
             'gpmB', 'apgM', 'serA', 'serC', 'serB', 'thrH', 'psp', 'porA',
             'porB', 'porD', 'porG', 'pps', 'mdh-K00024'),
  rump = c('hxlB', 'hxlA', 'pfkA', 'pfkB', 'pfkC', 'pfp', 'fbaA', 'fbaB',
           'glpX', 'fbp3', 'fae-hps', 'fbp-SEBP'),


  # nitrogen cycles
  nitrification = c('amoA_A', 'amoB_A', 'amoC_A', 'amoA_B', 'amoB_B', 'amoC_B',
                    'hao', 'nxrA', 'nxrB'),
  denitrification = c('napA', 'napB', 'napC', 'narG', 'narH', 'narJ', 'narI',
                      'nirK', 'nirS', 'norB', 'norC', 'nosZ', 'narZ', 'narY', 'narV', 'narW'),
  assnitred = c('nasA', 'nasB', 'nirA', 'NR', 'narB', 'narC'),
  dissnitred = c('napA', 'napB', 'napC', 'narG', 'narH', 'narJ', 'narI', 'narZ',
                 'narY', 'narV', 'narW', 'nirB', 'nirD', 'nrfA', 'nrfB', 'nrfC',
                 'nrfD'),
  nitfix = c('anfG', 'nifD', 'nifH', 'nifK', 'nifW'),
  annamox = c('hzo', 'hzsA', 'hzsB', 'hzsC', 'hdh'),
  odegsyn = c('ureA', 'ureB', 'ureC', 'nao', 'nmo', 'gdh_K00260', 'gdh_K00261',
              'gdh_K00262', 'gdh_K15371', 'gs_K00264', 'gs_K00265', 'gs_K00266',
              'gs_K00284', 'glsA', 'glnA', 'asnB', 'ansB'),
  nitrogen_other = c('hcp', 'pmoA', 'pmoB', 'pmoC'),



  # phosphorus cycles
  pyruvate = c('pps', 'ppdK', 'pyk', 'pckG', 'ppc', 'pckA'),
  pentose = c('gdh', 'gcd', 'gnl', 'gntK', 'gnd', 'rpiA', 'prsA', 'deoB'),
  phosphotransferase = c('ptsI', 'ptsH'),
  ox_phosphorylation = c('ppk', 'ppa'),
  phosph_met = c('pepM', 'pphA', 'ppd', 'phnX', 'fomC', 'phpC', 'mpnS', 'phnG',
                 'phnH', 'phnI', 'phnK', 'phnL', 'phnM', 'phnJ', 'phnP', 'phnN',
                 'phnPP', 'phny', 'phnA', 'phnW', 'phnO', 'pbfA', 'phnY', 'phnZ'),
  two_comp = c('phoU', 'phoR', 'phoB', 'phoP', 'SenX3', 'RegX3', 'pgtC', 'pgtB', 'pgtA'),
  transporters = c('pgtP', 'pstS', 'pstC', 'pstA', 'pstB', 'pit', 'htxB', 'ptxA',
                   'ptxB', 'ptxC', 'phnD_phosphite', 'phnD', 'phnE', 'phnC',
                   'ugpB', 'ugpA', 'ugpE', 'ugpC', 'phnS', 'phnV', 'phnU',
                   'phnT', 'glpT', 'aepX', 'aepV', 'aepW', 'aepP', 'aepS'),
  org_phos_hyd = c('opd', 'pafA', 'phoA', 'phoD', 'phoX', 'phoN', 'aphA',
                   'phoC', 'olpA', 'phy', 'appA', 'ugpQ', 'glpQ'),
  phos_other = c('htxA', 'ptxD', 'lysR', 'phnR', 'phnF', 'phoH'),
  purine = c('purF', 'purD', 'purN', 'purT', 'purL', 'purS', 'purQ', 'purM',
             'purK', 'purE', 'ADE2', 'purC', 'purB', 'purH', 'purP', 'purO',
             'guaB', 'guaA', 'gmk', 'ushA', 'ndk', 'spoT', 'ppx', 'purA', 'adk'),
  pyrimidine = c('pyrE', 'pyrF', 'ushA', 'cmk', 'pyrH', 'ndk', 'pyrG', 'rtpR',
                 'nrdD', 'nrdA', 'nrdE', 'nrdB', 'nrdF', 'nrdJ', 'dcd', 'dut',
                 'thyA', 'tmk'),


  #sulfur cycles
  assulred = c('sat', 'sir', 'cysC', 'cysN', 'cysD', 'cysH', 'cysI', 'cysJ',
               'cysN', 'cysC', 'cysQ', 'nrnA'),
  dsro = c('aprA', 'aprB', 'dsrA', 'dsrB', 'dsrC', 'dsrD', 'dsrN', 'dsrT',
           'dsrE', 'dsrF', 'dsrH', 'dsrL', 'dsrM', 'dsrK', 'dsrJ', 'dsrO',
           'dsrP', 'qmoA', 'qmoB', 'qmoC', 'rdsr', 'sat'),
  sulred = c('asrA', 'asrB', 'asrC', 'fsr', 'hydA', 'hydB', 'hydD', 'hydG',
             'mccA', 'otr', 'psrA', 'psrB', 'psrC', 'rdlA', 'shyA', 'shyB',
             'shyC', 'shyD', 'sreA', 'sreB', 'sreC', 'sudA', 'sudB', 'ttrA',
             'ttrB', 'ttrC'),
  SOX = c('soxA', 'soxX', 'soxB', 'soxC', 'soxD', 'soxY', 'soxZ'),
  sulfur_ox = c('doxA', 'doxD', 'fccA', 'fccB', 'glpE', 'soeA', 'soeB', 'soeC',
                'sorA', 'sorB', 'sqr', 'sseA', 'tsdA', 'tsdB'),
  sulfur_dis = c('phsA', 'phsB', 'phsC', 'tetH', 'sor'),
  org_sul_trans = c('acuI', 'acuN', 'acuK', 'betA', 'betB', 'betC', 'comA',
                    'comB', 'comC', 'comD', 'comE', 'dddA', 'dddC', 'dddD',
                    'dddK', 'dddL', 'dddP', 'dddQ', 'dddW', 'dddY', 'dddT',
                    'ddhA', 'ddhB', 'ddhC', 'dmdA', 'dmdB', 'dmdC', 'dmdD',
                    'dmoA', 'dmsA', 'dmsB', 'dmsC', 'dsyB', 'gdh', 'hpsN',
                    'hpsO', 'hpsP', 'iseJ', 'isfD', 'mddA', 'mdh', 'mtsA',
                    'mtsB', 'prpE', 'pta', 'sfnG', 'slcC', 'slcD', 'sqdB',
                    'sqdD', 'sqdX', 'tauX', 'tauY', 'tmm', 'toa', 'tpa', 'yihQ'),
  in_or_sul = c('cuyA', 'cysE', 'cysK', 'cysM', 'cysO', 'hdrA', 'hdrB', 'hdrC',
                'hdrD', 'hdrE', 'mccB', 'metA', 'metB', 'metC', 'metX', 'metY',
                'metZ', 'msmA', 'msmB', 'mtoX', 'ssuD', 'ssuE', 'suyA', 'suyB',
                'tauD', 'tbuB', 'tbuC', 'tmoC', 'tmoF', 'touC', 'touF', 'xsc'),
  sul_other = c('cuyZ', 'cysA', 'cysP', 'cysU', 'cysW', 'cysZ', 'hpsK', 'hpsL',
                'hpsM', 'iseK', 'iseL', 'iseM', 'sbp', 'sgpA', 'sgpB', 'sgpC',
                'soxL', 'ssuA', 'ssuB', 'ssuC', 'sulP', 'tauA', 'tauB', 'tauC',
                'tauE', 'tauZ', 'tusA', 'tusB', 'tusC', 'tusD', 'tusE')

)

short_pathways <- c(
  "cenmetpat", "hydro_met", "aceti_met", "methyl_met", "aom", "oxid_met_c1",
  "oxid_formaldehyde", "oxid_formate", "serine", "rump",

  # nitrogen cycles
  "nitrification", "denitrification", "assnitred", "dissnitred", "nitfix",
  "annamox", "odegsyn", "nitrogen_other",

  # phosphorus cycles
  "pyruvate", "pentose", "phosphotransferase", "ox_phosphorylation",
  "phosph_met", "two_comp", "transporters", "org_phos_hyd", "phos_other",
  "purine", "pyrimidine",

  # sulfur cycles
  "assulred", "dsro", "sulred", "SOX", "sulfur_ox", "sulfur_dis",
  "org_sul_trans", "in_or_sul", "sul_other"
)

long_pathways <- c(
  "Central methanogenic pathway", "Hydrogenotrophic methanogenesis",
  "Aceticlastic methanogenesis", "Methylotrophic methanogenesis",
  "Anaerobic oxidation of methane (AOM)", "Oxidation of methane and C1 compounds",
  "Oxidation of formaldehyde", "Oxidation of formate", "Serine cycle", "RuMP cycle",

  "Nitrification", "Denitrication", "Assimilatory nitrate reduction",
  "Dissimilatory nitrate reduction", "Nitrogen fixation", "Annamox",
  "Organic degradation and synthesis", "Related Nitrogen genes",

  "Pyruvate metabolism", "Pentose phosphate pathway", "Phosphotransferase system",
  "Oxidative phosphorylation", "Phosphonate and phosphinate metabolism",
  "Two-component system", "Transporters", "Organic phosphoester hydrolysis",
  "Related phosphorus genes", "Purine metabolism", "Pyrimidine metabolism",

  "Assimilatory sulphate reduction", "Dissimilatory sulphur reduction and oxidation",
  "Sulphur reduction", "SOX systems", "Sulphur oxidation", "Sulphur disproportionation",
  "Organic sulphur transformation",
  "Linkages between inorganic and organic sulphur transformation",
  "Related sulphur genes"
)

stopifnot(length(short_pathways) == length(long_pathways))
stopifnot(setequal(short_pathways, names(lop)))

pathwaymap <- data.frame(pathway = short_pathways, long_pathway = long_pathways, stringsAsFactors = FALSE)
