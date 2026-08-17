#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# KERRIA LACCA NOCICEPTION-RELATED ION CHANNEL PIPELINE
#
# Purpose:
#   Identify and validate the nociception-related ion-channel gene families analyzed by Goldberg et al. (2024):
#   pain, NompC, TRPM, Piezo, Pkd2, ppk/rpk/ppk26, and TRPA1.
#
# Reference assembly:
#   Kerria lacca GCA_045014175.1
#
# Main software:
#   eggNOG-mapper v2.1.15 + eggNOG v5.0.2
#   DIAMOND
#   NCBI BLAST+
#   MAFFT
#   IQ-TREE
#
# This script contains the final analytical path used for K. lacca.
################################################################################

# -----------------------------------------------------------------------------
# 0. USER SETTINGS
# -----------------------------------------------------------------------------

PROJECT="${HOME}/K_lacca_gene_families"
ASSEMBLY="GCA_045014175.1"

GENOME="${PROJECT}/data/K_lacca_GCA_045014175.1/ncbi_dataset/data/${ASSEMBLY}/${ASSEMBLY}_ASM4501417v1_genomic.fna"
GFF="${PROJECT}/data/K_lacca_GCA_045014175.1/ncbi_dataset/data/${ASSEMBLY}/genomic.gff"
CDS="${PROJECT}/data/K_lacca_GCA_045014175.1/ncbi_dataset/data/${ASSEMBLY}/cds_from_genomic.fna"
PROT="${PROJECT}/data/upload/K_lacca_GCA_045014175.1_proteins.faa"

EGGNOG_DB="${PROJECT}/data/eggnog_v5_db"
ANN="${PROJECT}/results/eggnog_v5_full/K_lacca_GCA_045014175_1.emapper.annotations"

mkdir -p \
  "${PROJECT}/results/eggnog_v5_full" \
  "${PROJECT}/results/piezo_validation" \
  "${PROJECT}/results/pkd2_validation" \
  "${PROJECT}/results/ppk_phylogeny" \
  "${PROJECT}/data/diamond" \
  "${PROJECT}/data/blast" \
  "${PROJECT}/data/references" \
  "${PROJECT}/logs"

cd "${PROJECT}"

# -----------------------------------------------------------------------------
# 1. SOFTWARE ENVIRONMENTS
# -----------------------------------------------------------------------------
# The original analysis used three Conda environments:
#   eggnog_v5   : eggNOG-mapper v2.1.15, Python 3.8
#   blast_tools : NCBI BLAST+
#   phylo_tools : MAFFT + IQ-TREE
#
# From a new shell:
# source "$HOME/miniforge3/etc/profile.d/conda.sh"
#
# Activate the environment indicated at the beginning of each section.
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# 2. INPUT CHECKS
# -----------------------------------------------------------------------------

echo "Protein count:"
grep -c '^>' "${PROT}"
echo "Expected: 10696 proteins"

# -----------------------------------------------------------------------------
# 3. EGGNOG-MAPPER ANNOTATION
# -----------------------------------------------------------------------------
# Environment: eggnog_v5
# conda activate eggnog_v5

emapper.py \
  -m diamond \
  -i "${PROT}" \
  --itype proteins \
  --data_dir "${EGGNOG_DB}" \
  -o K_lacca_GCA_045014175_1 \
  --output_dir results/eggnog_v5_full \
  --cpu 4 \
  2>&1 | tee logs/emapper_K_lacca_full.log

# -----------------------------------------------------------------------------
# 4. TARGET ORTHOGROUP SCREEN
# -----------------------------------------------------------------------------
# Goldberg et al. target orthogroups:
# pain            41TMZ
# NompC           41U79
# TRPM            41UGV
# Piezo           41UWY
# Pkd2            41VXR
# ppk/rpk/ppk26   41X6F
# TRPA1           41XIN

echo -e "Orthogroup\tCount" > results/K_lacca_target_orthogroup_counts.tsv

for og in 41TMZ 41U79 41UGV 41UWY 41VXR 41X6F 41XIN; do
  n=$(awk -F'\t' -v og="${og}" '
    $1 !~ /^#/ && index($5, og) {count++}
    END {print count+0}
  ' "${ANN}")
  echo -e "${og}\t${n}" | tee -a results/K_lacca_target_orthogroup_counts.tsv
done

# Expected counts:
# 41TMZ  1
# 41U79  1
# 41UGV  1
# 41UWY  2
# 41VXR  0
# 41X6F  2
# 41XIN  1

awk -F'\t' '
BEGIN {
  OFS="\t";
  print "Orthogroup","Query","Seed_ortholog","Evalue","Score","Preferred_name","PFAMs"
}
$1 !~ /^#/ {
  split("41TMZ 41U79 41UGV 41UWY 41VXR 41X6F 41XIN", ogs, " ")
  for (i in ogs) {
    if (index($5, ogs[i])) {
      print ogs[i],$1,$2,$3,$4,$9,$21
    }
  }
}' "${ANN}" > results/K_lacca_target_candidates.tsv

# -----------------------------------------------------------------------------
# 5. STRAIGHTFORWARD CANDIDATES
# -----------------------------------------------------------------------------
# Final candidates:
# pain   : XKL66442.1
# NompC  : XKL61404.1
# TRPM   : XKL64618.1
# TRPA1  : XKL60332.1

grep -E '^#query|^XKL66442\.1|^XKL61404\.1|^XKL64618\.1|^XKL60332\.1' \
  "${ANN}" | cut -f1,2,3,4,8,9,21 \
  > results/K_lacca_simple_candidates_annotation.tsv

# Independent TRPM copy-number check.
awk -F'\t' '$21 ~ /TRPM_tetra/ {
  print $1,$5,$9,$21
}' OFS='\t' "${ANN}" \
  > results/K_lacca_TRPM_tetra_hits.tsv

# Expected: only XKL64618.1

# -----------------------------------------------------------------------------
# 6. BUILD LOCAL K. LACCA DIAMOND DATABASE
# -----------------------------------------------------------------------------
# Environment: eggnog_v5
# conda activate eggnog_v5

DIAMOND="$(python -c 'import eggnogmapper, os; print(os.path.join(os.path.dirname(eggnogmapper.__file__), "bin", "diamond"))')"

"${DIAMOND}" makedb \
  --in "${PROT}" \
  -d data/diamond/K_lacca_proteome

# -----------------------------------------------------------------------------
# 7. PIEZO VALIDATION
# -----------------------------------------------------------------------------

awk '
/^>/ {
  keep = ($1==">XKL69283.1" || $1==">XKL69451.1")
}
keep
' "${PROT}" > results/piezo_candidates.faa

# Environment: blast_tools
# conda activate blast_tools
blastp \
  -query results/piezo_candidates.faa \
  -subject results/piezo_candidates.faa \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  > results/piezo_candidates_vs_each_other.tsv

awk '$1 != $2' results/piezo_candidates_vs_each_other.tsv \
  > results/piezo_candidates_pairwise_nonself.tsv

grep -n -E 'XKL69283\.1|XKL69451\.1' "${GFF}" \
  > results/piezo_validation/K_lacca_Piezo_GFF_hits.txt

grep -E '^#query|^XKL69283\.1|^XKL69451\.1' "${ANN}" \
  | cut -f1,21 \
  > results/piezo_validation/K_lacca_Piezo_PFAMs.tsv

# Environment: eggnog_v5
# conda activate eggnog_v5
"${DIAMOND}" blastp \
  --query results/piezo_candidates.faa \
  --db "${EGGNOG_DB}/eggnog_proteins.dmnd" \
  --out results/piezo_validation/Klacca_Piezo_vs_eggnog.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore full_sseq \
  --very-sensitive \
  --max-target-seqs 50 \
  --evalue 1e-5 \
  --threads 4

awk '$2=="7029.ACYPI50514-PA" {
  print ">"$2
  print $7
  exit
}' results/piezo_validation/Klacca_Piezo_vs_eggnog.tsv \
  > data/references/Apisum_Piezo_ACYPI50514-PA.faa

"${DIAMOND}" blastp \
  --query data/references/Apisum_Piezo_ACYPI50514-PA.faa \
  --db data/diamond/K_lacca_proteome \
  --out results/piezo_validation/Apisum_Piezo_vs_Klacca_proteome.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore \
  --very-sensitive \
  --max-target-seqs 20 \
  --evalue 1e-5 \
  --threads 4

# Key observed results:
# XKL69283.1  46.7% identity, 98.3% reference coverage, E=0
# XKL69451.1  24.9% identity, 88.0% reference coverage, E=2.49e-138
# Decision: retain both as independent Piezo-family loci.

# -----------------------------------------------------------------------------
# 8. PKD2 VALIDATION
# -----------------------------------------------------------------------------
# No protein was assigned to 41VXR.
# Reference FASTA files used:
#   data/references/Dmel_Pkd2_FBpp0079947.faa
#   data/references/Apisum_41VXR_ACYPI30693-PA.faa
#
# The D. melanogaster reference can be extracted from the FlyBase translation
# FASTA if it is available in data/references/.

if [[ -f data/references/dmel-all-translation-r6.68.fasta.gz ]]; then
  zcat data/references/dmel-all-translation-r6.68.fasta.gz | \
  awk '
  /^>/ {keep = ($1==">FBpp0079947")}
  keep
  ' > data/references/Dmel_Pkd2_FBpp0079947.faa
fi

# Environment: eggnog_v5
# conda activate eggnog_v5
"${DIAMOND}" blastp \
  --query data/references/Dmel_Pkd2_FBpp0079947.faa \
  --db data/diamond/K_lacca_proteome \
  --out results/pkd2_validation/Dmel_Pkd2_vs_K_lacca.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore \
  --very-sensitive \
  --max-target-seqs 20 \
  --evalue 1e-5 \
  --threads 4

"${DIAMOND}" blastp \
  --query data/references/Apisum_41VXR_ACYPI30693-PA.faa \
  --db data/diamond/K_lacca_proteome \
  --out results/pkd2_validation/Apisum_Pkd2_vs_K_lacca.tsv \
  --outfmt 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore \
  --very-sensitive \
  --max-target-seqs 20 \
  --evalue 1e-5 \
  --threads 4

# Environment: blast_tools
# conda activate blast_tools
makeblastdb \
  -in "${GENOME}" \
  -dbtype nucl \
  -parse_seqids \
  -out data/blast/K_lacca_genome

tblastn \
  -query data/references/Apisum_41VXR_ACYPI30693-PA.faa \
  -db data/blast/K_lacca_genome \
  -evalue 1e-3 \
  -max_target_seqs 100 \
  -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
  -out results/pkd2_validation/Apisum_Pkd2_vs_K_lacca_genome.tsv

tblastn \
  -query data/references/Apisum_41VXR_ACYPI30693-PA.faa \
  -db data/blast/K_lacca_genome \
  -matrix BLOSUM45 \
  -word_size 2 \
  -evalue 10 \
  -max_target_seqs 200 \
  -num_threads 4 \
  -outfmt "6 qseqid sseqid pident length qlen qstart qend sstart send evalue bitscore" \
  -out results/pkd2_validation/Apisum_Pkd2_vs_K_lacca_genome_sensitive.tsv

# Observed permissive hit:
# 60 aa, 25% identity, E=3.7, bitscore=33.7
# Decision: reject as background similarity.
# Final status: Pkd2 NOT DETECTED, not "confirmed absent".

# -----------------------------------------------------------------------------
# 9. PPK / RPK / PPK26 FAMILY
# -----------------------------------------------------------------------------
# eggNOG 41X6F candidates:
#   XKL60324.1
#   XKL60325.1
# Goldberg et al. treat ppk/rpk/ppk26 as one family.

awk '
/^>/ {
  keep = ($1==">XKL60324.1" || $1==">XKL60325.1")
}
keep
' "${PROT}" > results/ppk_family_candidates.faa

grep -n -E 'XKL60324\.1|XKL60325\.1' "${GFF}" \
  > results/ppk_phylogeny/K_lacca_41X6F_GFF_hits.txt

# Environment: blast_tools
# conda activate blast_tools
blastp \
  -query results/ppk_family_candidates.faa \
  -subject results/ppk_family_candidates.faa \
  -evalue 1e-5 \
  -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \
  > results/ppk_family_candidates_vs_each_other.tsv

awk '$1 != $2' results/ppk_family_candidates_vs_each_other.tsv \
  > results/ppk_phylogeny/K_lacca_41X6F_pairwise.tsv

# Identify all K. lacca ASC/DEG-ENaC proteins.
awk -F'\t' '$21 ~ /(^|,)ASC(,|$)/ {
  print $1,$5,$9
}' OFS='\t' "${ANN}" \
  > results/ppk_phylogeny/K_lacca_all_ASC.tsv

# Extract the 7 ASC proteins detected in the final analysis.
awk '
/^>/ {
  keep = ($1==">XKL60324.1" ||
          $1==">XKL60325.1" ||
          $1==">XKL62587.1" ||
          $1==">XKL62645.1" ||
          $1==">XKL63516.1" ||
          $1==">XKL65360.1" ||
          $1==">XKL69544.1")
}
keep
' "${PROT}" > results/ppk_phylogeny/Klacca_all_7_ASC.faa

# D. melanogaster references.
# Requires data/references/dmel-all-translation-r6.68.fasta.gz
if [[ -f data/references/dmel-all-translation-r6.68.fasta.gz ]]; then
  zcat data/references/dmel-all-translation-r6.68.fasta.gz | \
  awk '
  /^>/ {
    keep = ($1==">FBpp0080210" ||
            $1==">FBpp0293856" ||
            $1==">FBpp0305020")
  }
  keep
  ' > data/references/Dmel_ppk_rpk_ppk26.faa
fi

# FBpp0080210 -> ppk
# FBpp0293856 -> rpk
# FBpp0305020 -> ppk26

# A. pisum 41X6F references.
# Environment: eggnog_v5
# conda activate eggnog_v5
"${DIAMOND}" blastp \
  --query results/ppk_family_candidates.faa \
  --db "${EGGNOG_DB}/eggnog_proteins.dmnd" \
  --out results/ppk_family_vs_eggnog.tsv \
  --outfmt 6 qseqid sseqid pident length evalue bitscore full_sseq \
  --very-sensitive \
  --max-target-seqs 200 \
  --evalue 1e-5 \
  --threads 4

awk '
$2=="7029.ACYPI068872-PA" && !seen[$2]++ {print ">"$2; print $7}
$2=="7029.ACYPI29894-PA"  && !seen[$2]++ {print ">"$2; print $7}
$2=="7029.ACYPI35976-PA"  && !seen[$2]++ {print ">"$2; print $7}
' results/ppk_family_vs_eggnog.tsv \
  > data/references/Apisum_41X6F_three.faa

# Expanded dipteran reference set.
# ppk   -> 44XMS@7147
# rpk   -> 44X3M@7147
# ppk26 -> 44YWK@7147
python - <<'PY' > results/ppk_phylogeny/diptera_ppk_reference_ids.tsv
import sqlite3

con = sqlite3.connect("data/eggnog_v5_db/eggnog.db")
ogs = {
    "ppk":   "44XMS@7147",
    "rpk":   "44X3M@7147",
    "ppk26": "44YWK@7147"
}

for gene, og in ogs.items():
    for name, pname in con.execute(
        "SELECT name, pname FROM prots WHERE ogs LIKE ? ORDER BY name",
        (f"%{og}%",)
    ):
        print(gene, name, pname if pname else "-", sep="\t")

con.close()
PY

# Expected: 76 rows.

"${DIAMOND}" blastp \
  --query data/references/Dmel_ppk_rpk_ppk26.faa \
  --db "${EGGNOG_DB}/eggnog_proteins.dmnd" \
  --out results/ppk_phylogeny/Dmel3_vs_eggnog_all.tsv \
  --outfmt 6 qseqid sseqid evalue bitscore full_sseq \
  --very-sensitive \
  --max-target-seqs 0 \
  --evalue 1e-5 \
  --threads 4

awk '
NR==FNR {wanted[$2]=1; next}
($2 in wanted) && !seen[$2]++ {
  print ">"$2
  print $5
}
' \
  results/ppk_phylogeny/diptera_ppk_reference_ids.tsv \
  results/ppk_phylogeny/Dmel3_vs_eggnog_all.tsv \
  > results/ppk_phylogeny/diptera_ppk_76refs.faa

awk '
NR==FNR {group[$2]=$1; next}
/^>/ {
  id=substr($1,2)
  print ">" group[id] "_" id
  next
}
{print}
' \
  results/ppk_phylogeny/diptera_ppk_reference_ids.tsv \
  results/ppk_phylogeny/diptera_ppk_76refs.faa \
  > results/ppk_phylogeny/diptera_ppk_76refs_labeled.faa

cat \
  results/ppk_phylogeny/diptera_ppk_76refs_labeled.faa \
  data/references/Apisum_41X6F_three.faa \
  results/ppk_phylogeny/Klacca_all_7_ASC.faa \
  > results/ppk_phylogeny/ppk_extended_86seq.faa

awk '
/^>/ {
  id=substr($1,2)

  if (id=="7029.ACYPI29894-PA")
    print ">A_pisum_ACYPI29894"
  else if (id=="7029.ACYPI068872-PA")
    print ">A_pisum_ACYPI068872"
  else if (id=="7029.ACYPI35976-PA")
    print ">A_pisum_ACYPI35976"
  else if (id=="XKL60324.1")
    print ">K_lacca_XKL60324"
  else if (id=="XKL60325.1")
    print ">K_lacca_XKL60325"
  else if (id=="XKL62587.1")
    print ">K_lacca_XKL62587"
  else if (id=="XKL62645.1")
    print ">K_lacca_XKL62645"
  else if (id=="XKL63516.1")
    print ">K_lacca_XKL63516"
  else if (id=="XKL65360.1")
    print ">K_lacca_XKL65360_ppk13"
  else if (id=="XKL69544.1")
    print ">K_lacca_XKL69544"
  else
    print $0

  next
}
{print}
' results/ppk_phylogeny/ppk_extended_86seq.faa \
  > results/ppk_phylogeny/ppk_extended_86seq_labeled.faa

# Environment: phylo_tools
# conda activate phylo_tools
mafft \
  --localpair \
  --maxiterate 1000 \
  --thread 4 \
  results/ppk_phylogeny/ppk_extended_86seq_labeled.faa \
  > results/ppk_phylogeny/ppk_extended_86seq_aligned.faa

iqtree \
  -s results/ppk_phylogeny/ppk_extended_86seq_aligned.faa \
  -m MFP \
  -B 1000 \
  -alrt 1000 \
  -T 4 \
  --prefix results/ppk_phylogeny/ppk_extended_86seq_ml

# Key result:
# XKL60324.1 + XKL60325.1 form a strongly supported pair
# (SH-aLRT/UFBoot = 97.4/100), but individual ppk/rpk/ppk26 orthology
# was not confidently resolved.
# Decision: report 2 copies at the FAMILY level.

# -----------------------------------------------------------------------------
# 10. FINAL VALIDATED CANDIDATE FASTA
# -----------------------------------------------------------------------------

awk '
/^>/ {
  keep = ($1==">XKL60324.1" ||
          $1==">XKL60325.1" ||
          $1==">XKL60332.1" ||
          $1==">XKL61404.1" ||
          $1==">XKL64618.1" ||
          $1==">XKL66442.1" ||
          $1==">XKL69283.1" ||
          $1==">XKL69451.1")
}
keep
' "${PROT}" > results/K_lacca_nociception_candidates.faa

echo "Final candidate protein count:"
grep -c '^>' results/K_lacca_nociception_candidates.faa
# Expected: 8

# -----------------------------------------------------------------------------
# 11. FINAL SUMMARY TABLE
# -----------------------------------------------------------------------------

cat > results/K_lacca_nociception_summary.tsv <<'TABLE'
Family	Orthogroup	Candidate	Length_aa	Detected_copies	Status
pain	41TMZ	XKL66442.1	943	1	Detected
NompC	41U79	XKL61404.1	1796	1	Detected
TRPM	41UGV	XKL64618.1	1665	1	Detected
Piezo	41UWY	XKL69283.1	2496	2	Detected
Piezo	41UWY	XKL69451.1	1587	2	Detected_divergent_member
Pkd2	41VXR	NA	NA	0	Not_detected
ppk/rpk/ppk26	41X6F	XKL60324.1	449	2	Detected_family_level
ppk/rpk/ppk26	41X6F	XKL60325.1	553	2	Detected_family_level
TRPA1	41XIN	XKL60332.1	1203	1	Detected
TABLE

echo
echo "Pipeline complete."
echo "Final summary: results/K_lacca_nociception_summary.tsv"
echo "Final FASTA:   results/K_lacca_nociception_candidates.faa"
