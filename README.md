# antarctic_seep_incubations_peerJ

Analysis pipeline for Antarctic seep sediment incubations: 16S amplicon + shotgun metagenomics, methane/oxygen manipulation across depths and sites (Cinder Cones Seep, Jetty), 2022 and 2023 field seasons.

### Layout

```
scripts/                analysis scripts by step -- see scripts/README.md
updated_16s_run/         raw 16S processing, QIIME2 artifacts, ASV table, metadata
kaiju_output/            Kaiju taxonomy summaries + cached phyloseq
cycdb_diamond_output/    DIAMOND output vs cycDB + cached pathway counts
deseq_output/            DESeq2 pathway results + lollipop plots
figures/                 main-text figures
supplementary/           supplementary figures/tables
dic_rates.csv            DIC production rate data
```

`.gitignore` excludes raw reads (`*.fastq.gz`), raw DIAMOND tsvs, and QIIME2/SILVA reference files -- too large for git, regenerable from `scripts/`.

See `scripts/README.md` for what each script does.
