# Manuscript_Candida_Albicans

Correlating *Candida albicans* (and other gut taxa) with host mucosal
transcription in the HMP2 / IBDMDB cohort, using the
**Kraken Plus PFP (2024-10)** species-level reference.

This repository is the self-contained analysis pipeline behind the
Candida-vs-host-transcriptome results: it starts from a species CLR
table and a biopsy RNA-seq count matrix, and produces per-diagnosis
gene-wise correlations, HNF4A transcription-factor enrichment, GO
enrichment, BTNL3 / BTNL8 direct regressions, a week-window
sensitivity analysis, and two combined HNF4A + BTNL3 + BTNL8 taxon
rankings (one over a curated 13-species fungal panel, one over the
top-100 most abundant HD taxa plus a forced *C. albicans*).

> **TL;DR.** Clone, install the five Python dependencies listed below,
> then run `bash runme.sh`. Everything the pipeline needs — the
> species table, the host transcriptomics counts, the metadata, the
> Enrichr gene-set libraries, and a cached mygene.info annotation
> table — is already under `Data/`. No network access is required.

---

## Cohort and data

- **Cohort:** HMP2 / IBDMDB (SRA BioProject
  [PRJNA398089](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA398089)),
  longitudinal stool shotgun metagenomics + matched biopsy host
  transcriptomics from Crohn's disease (CD), ulcerative colitis (UC)
  and non-IBD control participants. Throughout the pipeline, "HD"
  (healthy donor) includes both the original HD and nonIBD labels.
- **Metagenomics:** species-level CLR table produced by running the
  HMP2 stool shotgun reads against the
  [Kraken 2 / Bracken "PFP" prebuilt index, 2024-10 release](https://benlangmead.github.io/aws-indexes/k2)
  from the Langmead lab. Zeros are treated as missing before CLR.
  The resulting `Run × species` table has 17,031 species.
- **Host transcriptomics:** biopsy RNA-seq count matrix from the HMP2
  data portal. Row labels are filtered through
  [mygene.info](https://mygene.info) (a BioThings service aggregating
  Ensembl / NCBI Entrez / UniProt / HGNC) to retain only rows that
  match a known human gene symbol or alias; ~33,052 genes survive
  and form the background for every enrichment test.
- **Gene-set libraries:** four GMT files from the Ma'ayan Lab
  [Enrichr](https://maayanlab.cloud/Enrichr) resource, bundled under
  `Data/`:
  - `ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt` — consensus
    ChIP-seq / ChIP-chip targets merged from ENCODE and ChEA; used
    for transcription-factor enrichment.
  - `GO_Biological_Process_2023.gmt`,
    `GO_Molecular_Function_2023.gmt`,
    `GO_Cellular_Component_2023.gmt` — 2023 release of the GO
    namespaces, used for negative-correlation enrichment.

---

## Repository layout

```
Manuscript_Candida_Albicans/
├── runme.sh                         # orchestrates the 8-step pipeline
├── methods.txt                      # Nature Microbiology-style methods
├── README.md                        # this file
│
├── Data/                            # all input data, fully self-contained
│   ├── species_table_clr_202410.tsv     # Kraken PFP 2024-10 CLR, Run × species
│   ├── PRJNA398089_MGX_merged.tsv       # metagenomics run metadata
│   ├── hmp2_metadata_2018-08-20.csv     # HMP2 portal metadata
│   ├── host_tx_counts.tsv               # raw biopsy RNA-seq counts
│   ├── host_tx_counts_filtered.tsv      # mygene.info-filtered counts (33,052 genes)
│   ├── gene_annotations.tsv             # cached mygene.info annotation table
│   ├── ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt
│   ├── GO_Biological_Process_2023.gmt
│   ├── GO_Cellular_Component_2023.gmt
│   └── GO_Molecular_Function_2023.gmt
│
├── _figutil.py                      # helper: save every plot as .png + .pdf
├── filter_known_genes.py            # Step 1 — mygene.info filter / cache
├── correlate_candida_genes.py       # Step 2 — Candida vs. gene Pearson, per diagnosis
├── enrichr_enrichment.py            # Steps 3–4 — TF + GO enrichment
├── btnl_analysis.py                 # Step 5 — BTNL3 / BTNL8 direct analysis
├── hnf4a_week_sensitivity.py        # Step 6 — HNF4A enrichment vs. ±week window
├── hnf4a_btnl_fungal_ranking.py     # Step 7 — 13-fungi combined ranking
├── hnf4a_btnl_taxa_ranking.py       # Step 8 — top-100 HD taxa combined ranking
│
├── Results/                         # per-diagnosis correlations + enrichment
│   ├── HD/{Tables,Figures}/
│   ├── CD/{Tables,Figures}/
│   ├── UC/{Tables,Figures}/
│   └── BTNL/                        # BTNL3 / BTNL8 scatters + week sensitivity
├── Results-HNF4A-sensitivity/       # HNF4A enrichment vs. ±week window (HD/CD/UC)
├── Results-HNF4A-fungi/              # 13 curated fungi — HNF4A + BTNL3 + BTNL8
│   └── methods.txt                  #   (including fungal curation rationale)
└── Results-HNF4A-taxa/               # top-100 HD taxa + Candida — same scoring
    └── methods.txt
```

Every figure the pipeline produces is written as **both** `.png`
(raster, 200 dpi) and `.pdf` (vector), via the shared helper
`_figutil.save_fig`.

---

## Requirements

- Python ≥ 3.9
- `pandas`, `numpy`, `scipy`, `matplotlib`, `tqdm`

Install with:

```bash
pip install pandas numpy scipy matplotlib tqdm
```

No other tooling (no R, no Snakemake, no Docker) is required.

### Network access

The pipeline runs **entirely offline** as shipped: `Data/gene_annotations.tsv`
is a cached mygene.info annotation table (see
[How the mygene.info filter works](#how-the-mygeneinfo-filter-works)
below), and the Enrichr libraries are bundled GMT files.

If you want to force a fresh mygene.info query — for example after
updating the raw `host_tx_counts.tsv` — delete
`Data/gene_annotations.tsv` *and* `Data/host_tx_counts_filtered.tsv`
and re-run, or pass `--refresh` to `filter_known_genes.py`.

---

## Running the pipeline

From the repository root:

```bash
bash runme.sh
```

Total runtime on a laptop is roughly **25–40 minutes**, dominated by
step 2 (per-gene `scipy.stats.linregress` for 33 k genes × 3
diagnoses) and step 6 (week-window sweep). Output is streamed to
stdout; redirect to a log file if you want to keep it:

```bash
bash runme.sh 2>&1 | tee runme.log
```

### What the 8 steps do

| Step | Script                              | What it produces |
|------|-------------------------------------|------------------|
| 1    | `filter_known_genes.py`             | `Data/host_tx_counts_filtered.tsv` (restricted to 33 k known human genes) and `Data/gene_annotations.tsv` (the mygene.info cache). Skipped if the filtered file already exists. |
| 2    | `correlate_candida_genes.py`        | Per-diagnosis Pearson `r`, slope and *p*-value of every host gene vs. *C. albicans* CLR, plus scatter plots for a curated gene panel (TLR / IL / NOD2 / HNF4A / HNF4G). Results under `Results/{HD,CD,UC}/`. |
| 3    | `enrichr_enrichment.py` (TF)        | Fisher's-exact TF enrichment of the ENCODE_and_ChEA_Consensus library against host genes *positively* correlated with *C. albicans* at `p < 0.05`, per diagnosis. |
| 4    | `enrichr_enrichment.py` (GO BP/MF/CC)| Same, for the three 2023 GO namespaces, on *negatively* correlated genes. |
| 5    | `btnl_analysis.py`                  | Direct Pearson regression of BTNL3 and BTNL8 vs. *C. albicans* CLR; scatter plots plus a week-window sensitivity sweep. Outputs under `Results/BTNL/`. |
| 6    | `hnf4a_week_sensitivity.py`         | HNF4A ENCODE enrichment *p*-value and BH-adjusted *p*-value as a function of the temporal matching window (±1 to ±52 weeks), for HD / CD / UC. Confirms the HNF4A signal is not an artefact of the default ±4-week window. |
| 7    | `hnf4a_btnl_fungal_ranking.py`      | Combined HNF4A + BTNL3 + BTNL8 score over **13 curated human-associated fungi**, both uncapped and capped at `p = 1e-10`. Bar plots + TSVs under `Results-HNF4A-fungi/`. |
| 8    | `hnf4a_btnl_taxa_ranking.py`        | Same combined score over the **top-100 most abundant HD species plus *C. albicans*** (101 taxa total), using a vectorised per-taxon Pearson so all 33 k genes × 101 taxa fit in RAM. Outputs under `Results-HNF4A-taxa/`. |

### Key parameters

All are constants at the top of `runme.sh` / each script:

| Parameter        | Value                                              |
|------------------|----------------------------------------------------|
| `NUM_WEEKS`      | `4` — temporal matching window `|w_mgx − w_tx| ≤ 4`|
| Diagnoses        | `HD CD UC` (nonIBD merged into HD)                 |
| Target species   | `S__Candida albicans`                              |
| Corr. *p* threshold for query genes | `0.05`                          |
| TF of interest   | `HNF4A ENCODE`                                     |
| BTNL genes       | `BTNL3`, `BTNL8`                                   |
| Top-N taxa panel | `100` + forced *C. albicans*                       |
| Capped *p* floor | `1e-10`                                            |

---

## The combined score in one paragraph

For each candidate taxon we record three *p*-values: `HNF4A_pval`
(one-sided Fisher's exact enrichment of the HNF4A ENCODE target set
among genes positively correlated with the taxon's CLR at
uncorrected `p < 0.05`), and `BTNL3_pval` / `BTNL8_pval` (two-sided
Pearson *t*-test *p*-values for the direct regression of the two
BTNL genes against the taxon's CLR). The three are combined in log
space:

```
combined_geomean = exp( mean_i { log( max(p_i, 1e-300) ) } )
```

Because Fisher's exact can return HNF4A *p*-values as small as
`1e-55` for abundant Firmicutes, a single extreme hit can dominate
the geometric mean. As a complementary ranking that rewards taxa
whose **three signals are *jointly* strong** rather than those with
a single enormous HNF4A hit, we additionally cap every input
*p*-value at a floor of `1e-10` before taking the geometric mean;
both "uncapped" and "capped" rankings are reported.

---

## How the mygene.info filter works

We use [mygene.info](https://mygene.info) purely as a **gene-name
sanity filter**: row labels in the raw `host_tx_counts.tsv` are a
mixture of current HGNC symbols, retired symbols, Ensembl
scaffold/clone IDs with version suffixes, pseudogenes and non-
standard aliases. `filter_known_genes.py` batches unique row labels
1000 at a time into `https://mygene.info/v3/query` with
`scopes=["symbol","alias"]` and `species=human`, records whether
each label maps to a known human gene (and its `type_of_gene`), and
drops rows that don't match anything. The surviving ~33,052 rows
form the background for every downstream enrichment test.

The full annotation response is cached at `Data/gene_annotations.tsv`
so the pipeline can run offline thereafter. Passing `--refresh`
forces a re-query.

---

## Reproducibility notes

- The exact Kraken reference used was the Langmead lab prebuilt
  `k2_pluspfp_2024-10` database. Re-running against a different
  Kraken index will change the species table and therefore
  everything downstream.
- `scipy.stats.linregress` (correlate_candida_genes) and the
  vectorised numpy-matmul Pearson (hnf4a_btnl_taxa_ranking) are
  numerically equivalent to at least 10 significant figures on this
  data, and we have verified agreement for common taxa.
- Random state is irrelevant: no step uses stochastic sampling,
  bootstrapping or random permutation.
- Per-folder `methods.txt` files in `Results-HNF4A-fungi/` and
  `Results-HNF4A-taxa/` document the species-panel construction and
  output schema for those modules specifically. The root
  `methods.txt` is written in Nature Microbiology house style and
  covers the full pipeline.

---

## Citation

If you use this pipeline or its outputs, please cite the upstream
resources it depends on:

- **HMP2 / IBDMDB:** Lloyd-Price J. *et al.*, *Nature* **569**, 655–662
  (2019). BioProject PRJNA398089.
- **Kraken 2 / Bracken PFP index:** Wood D.E. *et al.*, *Genome Biology*
  **20**, 257 (2019); Lu J. *et al.*, *PeerJ Computer Science* **3**,
  e104 (2017); prebuilt index from the Langmead lab
  (https://benlangmead.github.io/aws-indexes/k2).
- **mygene.info / BioThings:** Xin J. *et al.*, *Genome Biology* **17**,
  91 (2016); Wu C. *et al.*, *Nucleic Acids Research* **41**, D561
  (2013).
- **Enrichr (Ma'ayan Lab):** Chen E.Y. *et al.*, *BMC Bioinformatics*
  **14**, 128 (2013); Kuleshov M.V. *et al.*, *Nucleic Acids Research*
  **44**, W90 (2016); Xie Z. *et al.*, *Current Protocols* **1**, e90
  (2021).
- **ENCODE/ChEA consensus ChIP-X:** as distributed by Enrichr.
- **Gene Ontology:** The Gene Ontology Consortium, *Nucleic Acids
  Research* **51**, D1257 (2023).

---

## License

Code in this repository is released under the MIT License. The
bundled Enrichr GMT libraries and the HMP2 metadata retain their
original upstream licensing terms; please refer to the Enrichr and
HMP2 data-use policies before redistributing those files.
