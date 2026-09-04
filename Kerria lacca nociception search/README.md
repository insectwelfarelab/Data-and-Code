# Kerria lacca nociception-related ion channel pipeline

This repository documents the workflow used to identify nociception-related ion channel gene families in the *Kerria lacca* genome assembly **GCA_045014175.1**.

The original analysis follows the comparative gene-family framework of Goldberg et al. (2024). An expanded analysis additionally screened three nociception-related targets: Straightjacket (stj), TRPA5, and Pyrexia (Pyx).

## Target families

| Family        | eggNOG v5 orthogroup |
|---------------|----------------------|
| pain          | 41TMZ                |
| NompC         | 41U79                |
| TRPM          | 41UGV                |
| Piezo         | 41UWY                |
| Pkd2          | 41VXR                |
| ppk/rpk/ppk26 | 41X6F                |
| TRPA1         | 41XIN                |

### Expanded targets

- Straightjacket (stj)
- TRPA5
- Pyrexia (Pyx)

## Required software

-   eggNOG-mapper v2.1.15
-   eggNOG v5.0.2 database
-   DIAMOND
-   NCBI BLAST+
-   MAFFT v7
-   IQ-TREE
-   Python 3 with `sqlite3`

The original analysis used three Conda environments:

-   `eggnog_v5`
-   `blast_tools`
-   `phylo_tools`

## Required inputs

The main script expects:

-   *K. lacca* genome assembly GCA_045014175.1
-   corresponding GFF3 annotation
-   predicted protein FASTA
-   eggNOG v5.0.2 database
-   FlyBase *D. melanogaster* translation FASTA for selected reference proteins
-   `Apisum_41VXR_ACYPI30693-PA.faa` for Pkd2 validation

The expanded analysis additionally uses curated insect reference proteins for Straightjacket and TRPA-family comparisons. For Pyrexia, 15 Pyx and 15 taxonomically matched Waterwitch (Wtrw) proteins were curated from Table S1 of Liénard et al. (2024).

Large genome and database files should not be committed to GitHub.

## Run

Original seven-family analysis:

``` bash
chmod +x scripts/K_lacca_nociception_pipeline.sh
bash scripts/K_lacca_nociception_pipeline.sh
```

Expanded analysis:

``` bash
chmod +x scripts/K_lacca_expanded_analysis_pipeline.sh
bash scripts/K_lacca_expanded_analysis_pipeline.sh
```

The original script is organized step-by-step and contains the final commands used for:

1.  eggNOG annotation
2.  target orthogroup screening
3.  domain validation
4.  Piezo validation
5.  Pkd2 protein- and genome-level validation
6.  ppk/rpk/ppk26 phylogenetic validation
7.  final candidate export

The expanded analysis uses eggNOG/PFAM annotation evidence, BLASTP searches against curated insect reference proteins, MAFFT alignments, and maximum-likelihood phylogenetic analyses with IQ-TREE. For TRPA5 and Pyrexia, genome-level TBLASTN searches were also performed to test for potentially unannotated loci.

## Final result

| Family        | Candidate(s)           | Detected copies |
|---------------|------------------------|----------------:|
| pain          | XKL66442.1             |               1 |
| NompC         | XKL61404.1             |               1 |
| TRPM          | XKL64618.1             |               1 |
| Piezo         | XKL69283.1; XKL69451.1 |               2 |
| Pkd2          | not detected           |               0 |
| ppk/rpk/ppk26 | XKL60324.1; XKL60325.1 |               2 |
| TRPA1         | XKL60332.1             |               1 |

Expanded analysis:

- Straightjacket: 1
- TRPA5: 0 detected
- Pyrexia: 0 detected

`Pkd2`, `TRPA5`, and `Pyrexia` are reported as **not detected** rather than definitively absent from the species.

The two `41X6F` proteins are reported at the **ppk/rpk/ppk26 family level** because individual orthology to ppk, rpk, or ppk26 was not confidently resolved.

## References

Goldberg JK, Godfrey RK, Barrett M. 2024. A long-read draft assembly of the Chinese mantis (Mantodea: Mantidae: *Tenodera sinensis*) genome reveals patterns of ion channel gain and loss across Arthropoda. *G3: Genes\|Genomes\|Genetics* 14:jkae062. <https://doi.org/10.1093/g3journal/jkae062>

Liénard MA, Baez-Nieto D, Tsai C-C, Valencia-Montoya WA, Werin B, Johanson U, Lassance J-M, Pan JQ, Yu N, Pierce NE. 2024. TRPA5 encodes a thermosensitive ankyrin ion channel receptor in a triatomine insect. *iScience* 27:109541. <https://doi.org/10.1016/j.isci.2024.109541>
