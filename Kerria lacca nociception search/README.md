This pipeline is the reference document for the work described here: https://forum.effectivealtruism.org/posts/2JdY8vuMAYcbqX5p9/research-report-a-genetic-basis-for-potential-nociception-in?sharePopup=true

# Kerria lacca nociception-related ion channel pipeline

This repository documents the workflow used to identify nociception-related ion channel gene families in the *Kerria lacca* genome assembly **GCA_045014175.1**.

The analysis follows the comparative gene-family framework of Goldberg et al. (2024).

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

Large genome and database files should not be committed to GitHub.

## Run

``` bash
chmod +x scripts/K_lacca_nociception_pipeline.sh
bash scripts/K_lacca_nociception_pipeline.sh
```

The script is organized by set-by-step and contains the final commands used for:

1.  eggNOG annotation
2.  target orthogroup screening
3.  domain validation
4.  Piezo validation
5.  Pkd2 protein- and genome-level validation
6.  ppk/rpk/ppk26 phylogenetic validation
7.  final candidate export

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

`Pkd2` is reported as **not detected**.

The two `41X6F` proteins are reported at the **ppk/rpk/ppk26 family level** because individual orthology to ppk, rpk, or ppk26 was not confidently resolved.

## Reference

Goldberg JK, Godfrey RK, Barrett M. 2024. A long-read draft assembly of the Chinese mantis (Mantodea: Mantidae: *Tenodera sinensis*) genome reveals patterns of ion channel gain and loss across Arthropoda. *G3: Genes\|Genomes\|Genetics* 14:jkae062. <https://doi.org/10.1093/g3journal/jkae062>
