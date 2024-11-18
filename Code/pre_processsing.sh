#!/bin/bash

# Function to print a progress message
progress() {
    local message="$1"
    local current_step="$2"
    local total_steps="$3"
    local width=50  # Width of the progress bar
    local percentage=$((current_step * 100 / total_steps))
    local num_completed=$((percentage * width / 100))
    local num_remaining=$((width - num_completed))
    local progress_bar="["
    progress_bar+="$(printf '=%.0s' $(seq 1 $num_completed))"
    progress_bar+="$(printf ' ' %.0s $(seq 1 $num_remaining))"
    progress_bar+="] $percentage%"
    echo -ne "\r$message\n$progress_bar"
}

# Total number of major steps in the script
total_steps=6  # Adjust this based on your script
current_step=1

# Step 1: Extract positions of the forward and reverse strand from the reference genome

progress "Step 1: Extracting positions from the reference genome" "$current_step" "$total_steps"
zcat Sus_scrofa.Sscrofa11.1.102.chr.gtf.gz | awk '($3=="gene")' | awk '($7=="+")'| awk 'BEGIN { OFS="\t" } {print$1, $4-200, $4+50, $10, $14}' | sed 's/"//g' | sed 's/;//g' > SS_11.1_v102_plus_strand.bed
zcat Sus_scrofa.Sscrofa11.1.102.chr.gtf.gz | awk '($3=="gene")' | awk '($7=="-")'| awk 'BEGIN { OFS="\t" } {print$1, $5-200, $5+50, $10, $14}' | sed 's/"//g' | sed 's/;//g' > SS_11.1_v102_minus_strand.bed
((current_step++))

# Step 2: Merge forward and reverse strand files

progress "Step 2: Merging forward and reverse strand files" "$current_step" "$total_steps"
cat SS_11.1_v102_plus_strand.bed SS_11.1_v102_minus_strand.bed > SS_11.1_v102.bed
((current_step++))

# Step 3: Join gene position file and DE file

progress "Step 3: Joining gene position file and DE file" "$current_step" "$total_steps"
join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.liverstage1.liverstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.liverstage2.liverstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.skinstage1.skinstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.skinstage2.skinstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.kidneystage1.kidneystage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.kidneystage2.kidneystage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.musclestage1.musclestage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.musclestage2.musclestage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.ileumstage1.ileumstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.ileumstage2.ileumstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.lungstage1.lungstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.lungstage2.lungstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.bed
((current_step++))

# Step 4: Intersect files (set - coordinates to 0 )

progress "Step 4: Intersecting files (setting - coordinates to 0 )" "$current_step" "$total_steps"
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.bed > GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.bed > GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed > GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed > GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.bed > GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.bed > GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.bed > GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.bed > GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.bed > GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.bed > GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.fixed.bed 
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.bed > GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.bed > GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.bed > GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.bed > GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.fixed.bed      

bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_liver_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_cerebellum_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_ileum_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_skin_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_lung_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.fixed.bed -b GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_muscle_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.fixed.bed -b GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_kidney_1_2_3.bed
((current_step++))


# Step 5: Filetr CpGs

progress "Step 5: Filetring CpGs" "$current_step" "$total_steps"
## Liver 30dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_2_CCAACACT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_3_CTCCAATC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_3.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_6_CCAACACT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_6.bed

## Liver 70dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_1_CCTATACC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_2_TAACGTCG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_3_CCTATACC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_3.bed

## Liver NB

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_M_1_CCACAACA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_M_2_GAAGACTG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_F_3_CCACAACA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_3.bed

## Hindbrain 30dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_2_CTAAGACC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_3_AATGACGC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_3.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_6_CTAAGACC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_6.bed

## Cerebellum 70dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_1_CAAGCCAA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_2_CCTGTCAA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_4_CAAGCCAA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_4.bed

## Cerebellum NB

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_M_1_GTGATCCA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_M_2_TCCAACTG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_F_3_AAGGACCA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_3.bed

## Skin 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_1_AACCGAAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_3_GTACCACA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_6_AACCGAAC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_6.bed

## Skin 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_1_TCACTCGA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_2_TCCTGGTA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_3_AGGTTCCT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_3.bed

## Skin NB

sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_F_3_CTATCCAC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_3.bed
sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_M_1_CTATCCAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_M_2_ACTGCACT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_2.bed

## Lung 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_2_CTCGAACA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_3_GTCTCATC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_6_CTCGAACA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_6.bed

## Lung 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_1_TCCATTGC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_2_AACACTGG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_3_TCCATTGC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_3.bed

## Lung NB

sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_F_4_ATAACGCC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_M_1_ATAACGCC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_M_2_TCTAGTCC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_2.bed

## Muscle 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_2_CCTTAGGT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_3_ACGATCAG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_6_CCTTAGGT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_6.bed

## Muscle 70dpf

sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_1_AAGTCCTC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_2_AGCAGACA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_4_AAGTCCTC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_4.bed

## Muscle NB

sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_F_4_ACGCTTCT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_M_1_ACGCTTCT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_M_2_GACATCTC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_2.bed

## Kidney 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_1_ACGGACTT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_3_GCCAGAAT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_6_ACGGACTT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_6.bed

## Kidney 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_1_CTGACTAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_2_CATCAACC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_3_GAACCTTC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_3.bed

## Kidney NB

sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_F_3_CCGGAATA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_3.bed
sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_M_1_CCGGAATA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_M_2_CTCGACTT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_2.bed

## Intestine 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_1_TCTAGGAG_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_3_ACTCTCCA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_6_TCTAGGAG_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_6.bed

## Intestine 70dpf

sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_1_AACGCCTT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_2_CGCAACTA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_4_AACGCCTT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_4.bed

## Intestine NB

sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_F_4_CCAAGTAG_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_M_1_CCAAGTAG_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_M_2_CTAGCTCA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_2.bed
((current_step++))

# Step 6: Finalize formatting

progress "Step 6: Finalize formatting" "$current_step" "$total_steps"
echo -e "chromosome\gene\tpval_1-2\tadjpval_1-2\tlogFC_1-2\tlogCPM_1-2\tpval_2-3\tadjpval_2-3\tlogFC_2-3\tlogCPM_2-3\tmeth_1\tmeth_2\tmeth_3" > header.txt

## Liver 30dpf

bedtools intersect -wo -a CpG_deg_liver_30dpf_POOL_2.bed -b CpG_deg_liver_30dpf_POOL_3.bed > ISEC_liver_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_liver_30dpf_2_3_temp.bed -b  CpG_deg_liver_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_liver_30dpf.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_POOLS_merged_liver_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_liver_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_liver_30dpf.bed > final_liver_30dpf.bed

## Liver 70dpf

bedtools intersect -wo -a CpG_deg_liver_70dpf_1.bed -b CpG_deg_liver_70dpf_2.bed > ISEC_liver_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_liver_70dpf_1_2_temp.bed -b  CpG_deg_liver_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_liver_70dpf.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_merged_liver_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_liver_70dpf.bed
cat header.txt temp_CpG_genes_merged_liver_70dpf.bed > final_liver_70dpf.bed

## Liver NB

bedtools intersect -wo -a CpG_deg_liver_NB_1.bed -b CpG_deg_liver_NB_2.bed > ISEC_liver_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_liver_NB_1_2_temp.bed -b  CpG_deg_liver_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_liver_NB.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_merged_liver_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_liver_NB.bed
cat header.txt temp_CpG_genes_merged_liver_NB.bed > final_liver_NB.bed

## Hindbrain 30dpf

bedtools intersect -wo -a CpG_deg_hindbrain_30dpf_POOL_2.bed -b CpG_deg_hindbrain_30dpf_POOL_3.bed > ISEC_hindbrain_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_hindbrain_30dpf_2_3_temp.bed -b  CpG_deg_hindbrain_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_hindbrain_30dpf.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_POOLS_merged_hindbrain_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_hindbrain_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_hindbrain_30dpf.bed > final_hindbrain_30dpf.bed

## Hindbrain 70dpf

bedtools intersect -wo -a CpG_deg_cerebellum_70dpf_1.bed -b CpG_deg_cerebellum_70dpf_2.bed > ISEC_cerebellum_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_cerebellum_70dpf_1_2_temp.bed -b  CpG_deg_cerebellum_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_cerebellum_70dpf.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_merged_cerebellum_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_cerebellum_70dpf.bed
cat header.txt temp_CpG_genes_merged_cerebellum_70dpf.bed > final_hindbrain_70dpf.bed

## Hindbrain NB

bedtools intersect -wo -a CpG_deg_cerebellum_NB_1.bed -b CpG_deg_cerebellum_NB_2.bed > ISEC_cerebellum_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_cerebellum_NB_1_2_temp.bed -b  CpG_deg_cerebellum_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_cerebellum_NB.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_merged_cerebellum_NB.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_cerebellum_NB.bed
cat header.txt temp_CpG_genes_merged_cerebellum_NB.bed > final_hindbrain_NB.bed

## Muscle 30dpf

bedtools intersect -wo -a CpG_deg_muscle_30dpf_POOL_2.bed -b CpG_deg_muscle_30dpf_POOL_3.bed > ISEC_muscle_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_muscle_30dpf_2_3_temp.bed -b  CpG_deg_muscle_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_muscle_30dpf.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_POOLS_merged_muscle_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_muscle_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_muscle_30dpf.bed > final_muscle_30dpf.bed

## Muscle 70dpf

bedtools intersect -wo -a CpG_deg_muscle_70dpf_1.bed -b CpG_deg_muscle_70dpf_2.bed > ISEC_muscle_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_muscle_70dpf_1_2_temp.bed -b  CpG_deg_muscle_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_muscle_70dpf.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_merged_muscle_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_muscle_70dpf.bed
cat header.txt temp_CpG_genes_merged_muscle_70dpf.bed > final_muscle_70dpf.bed

## Muscle NB

bedtools intersect -wo -a CpG_deg_muscle_NB_1.bed -b CpG_deg_muscle_NB_2.bed > ISEC_muscle_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_muscle_NB_1_2_temp.bed -b  CpG_deg_muscle_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_muscle_NB.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_merged_muscle_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_muscle_NB.bed
cat header.txt temp_CpG_genes_merged_muscle_NB.bed > final_muscle_NB.bed

## Kidney 30dpf

bedtools intersect -wo -a CpG_deg_kidney_30dpf_POOL_1.bed -b CpG_deg_kidney_30dpf_POOL_3.bed > ISEC_kidney_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_kidney_30dpf_1_3_temp.bed -b  CpG_deg_kidney_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_kidney_30dpf.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_POOLS_merged_kidney_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_kidney_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_kidney_30dpf.bed > final_kidney_30dpf.bed

## Kidney 70dpf

bedtools intersect -wo -a CpG_deg_kidney_70dpf_1.bed -b CpG_deg_kidney_70dpf_2.bed > ISEC_kidney_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_kidney_70dpf_1_2_temp.bed -b  CpG_deg_kidney_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_kidney_70dpf.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_merged_kidney_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_kidney_70dpf.bed
cat header.txt temp_CpG_genes_merged_kidney_70dpf.bed > final_kidney_70dpf.bed

## Kidney NB

bedtools intersect -wo -a CpG_deg_kidney_NB_1.bed -b CpG_deg_kidney_NB_2.bed > ISEC_kidney_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_kidney_NB_1_2_temp.bed -b  CpG_deg_kidney_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_kidney_NB.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_merged_kidney_NB.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_kidney_NB.bed
cat header.txt temp_CpG_genes_merged_kidney_NB.bed > final_kidney_NB.bed

## Lung 30dpf

bedtools intersect -wo -a CpG_deg_lung_30dpf_POOL_2.bed -b CpG_deg_lung_30dpf_POOL_3.bed > ISEC_lung_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_lung_30dpf_2_3_temp.bed -b  CpG_deg_lung_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_lung_30dpf.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_POOLS_merged_lung_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_lung_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_lung_30dpf.bed > final_lung_30dpf.bed

## Lung 70dpf

bedtools intersect -wo -a CpG_deg_lung_70dpf_1.bed -b CpG_deg_lung_70dpf_2.bed > ISEC_lung_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_lung_70dpf_1_2_temp.bed -b  CpG_deg_lung_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_lung_70dpf.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_merged_lung_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_lung_70dpf.bed
cat header.txt temp_CpG_genes_merged_lung_70dpf.bed > final_lung_70dpf.bed

## Lung NB

bedtools intersect -wo -a CpG_deg_skin_NB_1.bed -b CpG_deg_skin_NB_2.bed > ISEC_lung_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_lung_NB_1_2_temp.bed -b  CpG_deg_skin_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_lung_NB.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_merged_lung_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_lung_NB.bed
cat header.txt temp_CpG_genes_merged_lung_NB.bed > final_lung_NB.bed

## Skin 30dpf

bedtools intersect -wo -a CpG_deg_skin_30dpf_POOL_1.bed -b CpG_deg_skin_30dpf_POOL_3.bed > ISEC_skin_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_skin_30dpf_1_3_temp.bed -b  CpG_deg_skin_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_skin_30dpf.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_POOLS_merged_skin_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_skin_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_skin_30dpf.bed > final_skin_30dpf.bed

## Skin 70dpf

bedtools intersect -wo -a CpG_deg_skin_70dpf_1.bed -b CpG_deg_skin_70dpf_2.bed > ISEC_skin_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_skin_70dpf_1_2_temp.bed -b  CpG_deg_skin_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_skin_70dpf.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_merged_skin_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_skin_70dpf.bed
cat header.txt temp_CpG_genes_merged_skin_70dpf.bed > final_skin_70dpf.bed

## Skin NB

bedtools intersect -wo -a CpG_deg_skin_NB_1.bed -b CpG_deg_skin_NB_2.bed > ISEC_skin_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_skin_NB_1_2_temp.bed -b  CpG_deg_skin_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_skin_NB.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_merged_skin_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_skin_NB.bed
cat header.txt temp_CpG_genes_merged_skin_NB.bed > final_skin_NB.bed

## Intestine 30dpf

bedtools intersect -wo -a CpG_deg_intestine_30dpf_POOL_1.bed -b CpG_deg_intestine_30dpf_POOL_3.bed > ISEC_intestine_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_intestine_30dpf_1_3_temp.bed -b  CpG_deg_intestine_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_intestine_30dpf.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_POOLS_merged_intestine_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_intestine_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_intestine_30dpf.bed > final_intestine_30dpf.bed

## Intestine 70dpf

bedtools intersect -wo -a CpG_deg_intestine_70dpf_1.bed -b CpG_deg_intestine_70dpf_2.bed > ISEC_intestine_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_intestine_70dpf_1_2_temp.bed -b  CpG_deg_intestine_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_intestine_70dpf.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_merged_intestine_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_intestine_70dpf.bed
cat header.txt temp_CpG_genes_merged_intestine_70dpf.bed > final_intestine_70dpf.bed

## Intestine NB

bedtools intersect -wo -a CpG_deg_intestine_NB_1.bed -b CpG_deg_intestine_NB_2.bed > ISEC_intestine_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_intestine_NB_1_2_temp.bed -b  CpG_deg_intestine_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_intestine_NB.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_merged_intestine_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_intestine_NB.bed
cat header.txt temp_CpG_genes_merged_intestine_NB.bed > final_intestine_NB.bed

# Finally, update the progress bar to 100% when the script is done.
progress "Progress: Script completed" "$total_steps" "$total_steps"
exit 0n to print a progress message
progress() {
    local message="$1"
    local current_step="$2"
    local total_steps="$3"
    local width=50  # Width of the progress bar
    local percentage=$((current_step * 100 / total_steps))
    local num_completed=$((percentage * width / 100))
    local num_remaining=$((width - num_completed))
    local progress_bar="["
    progress_bar+="$(printf '=%.0s' $(seq 1 $num_completed))"
    progress_bar+="$(printf ' ' %.0s $(seq 1 $num_remaining))"
    progress_bar+="] $percentage%"
    echo -ne "\r$message\n$progress_bar"
}

# Total number of major steps in the script
total_steps=6  # Adjust this based on your script
current_step=1

# Step 1: Extract positions of the forward and reverse strand from the reference genome

progress "Step 1: Extracting positions from the reference genome" "$current_step" "$total_steps"
zcat Sus_scrofa.Sscrofa11.1.102.chr.gtf.gz | awk '($3=="gene")' | awk '($7=="+")'| awk 'BEGIN { OFS="\t" } {print$1, $4-200, $4+50, $10, $14}' | sed 's/"//g' | sed 's/;//g' > SS_11.1_v102_plus_strand.bed
zcat Sus_scrofa.Sscrofa11.1.102.chr.gtf.gz | awk '($3=="gene")' | awk '($7=="-")'| awk 'BEGIN { OFS="\t" } {print$1, $5-200, $5+50, $10, $14}' | sed 's/"//g' | sed 's/;//g' > SS_11.1_v102_minus_strand.bed
((current_step++))

# Step 2: Merge forward and reverse strand files

progress "Step 2: Merging forward and reverse strand files" "$current_step" "$total_steps"
cat SS_11.1_v102_plus_strand.bed SS_11.1_v102_minus_strand.bed > SS_11.1_v102.bed
((current_step++))

# Step 3: Join gene position file and DE file

progress "Step 3: Joining gene position file and DE file" "$current_step" "$total_steps"
join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.liverstage1.liverstage2.bed) | awk '{OFS="\t"} {print $1, $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.liverstage2.liverstage3.bed) | awk '{OFS="\t"} {print $1, $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed) | awk '{OFS="\t"} {print $1, $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.skinstage1.skinstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.skinstage2.skinstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.kidneystage1.kidneystage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.kidneystage2.kidneystage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.musclestage1.musclestage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.musclestage2.musclestage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.ileumstage1.ileumstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.ileumstage2.ileumstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.lungstage1.lungstage2.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.bed

join -1 4 -2 1 <(sort -k4 SS_11.1_v102.bed) <(sort -k1 diffres.pvals.diff.lungstage2.lungstage3.bed) | awk '{OFS="\t"} {print $2,$3,$4,$1,$5,$6,$7,$8,$9} ' > GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.bed
((current_step++))

# Step 4: Intersect files (set - coordinates to 0 )

progress "Step 4: Intersecting files (setting - coordinates to 0 )" "$current_step" "$total_steps"
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.bed > GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.bed > GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.bed > GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.bed > GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.bed > GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.bed > GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.bed > GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.bed > GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.bed > GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.bed > GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.fixed.bed 
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.bed > GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.bed > GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.bed > GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.fixed.bed
awk 'BEGIN {OFS="\t"} $2 < 0 { $2 = 0 } {print}' GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.bed > GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.fixed.bed      

bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.liverstage1.liverstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.liverstage2.liverstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_liver_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.cerebellumstage1.cerebellumstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.cerebellumstage2.cerebellumstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_cerebellum_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.ileumstage1.ileumstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.ileumstage2.ileumstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_ileum_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.skinstage1.skinstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.skinstage2.skinstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_skin_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.lungstage1.lungstage2.fixed.bed -b GENE_POS_diffres.pvals.diff.lungstage2.lungstage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_lung_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.musclestage1.musclestage2.fixed.bed -b GENE_POS_diffres.pvals.diff.musclestage2.musclestage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_muscle_1_2_3.bed
bedtools intersect -wo -a GENE_POS_diffres.pvals.diff.kidneystage1.kidneystage2.fixed.bed -b GENE_POS_diffres.pvals.diff.kidneystage2.kidneystage3.fixed.bed | awk '{OFS="\t"} {print$1,$2,$3,$4,$6,$7,$8,$9,$15,$16,$17,$18} ' > ISEC_deg_kidney_1_2_3.bed
((current_step++))


# Step 5: Filetr CpGs

progress "Step 5: Filetring CpGs" "$current_step" "$total_steps"
## Liver 30dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_2_CCAACACT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_3_CTCCAATC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_3.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_FT_30dpf_POOL_6_CCAACACT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_30dpf_POOL_6.bed

## Liver 70dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_1_CCTATACC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_2_TAACGTCG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_liver_FT_70dpf_3_CCTATACC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_70dpf_3.bed

## Liver NB

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_M_1_CCACAACA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_M_2_GAAGACTG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Liver_NB_F_3_CCACAACA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_liver_NB_3.bed

## Hindbrain 30dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_2_CTAAGACC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_3_AATGACGC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_3.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Hindbrain_FT_30dpf_POOL_6_CTAAGACC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_hindbrain_30dpf_POOL_6.bed

## Cerebellum 70dpf

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_1_CAAGCCAA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_2_CCTGTCAA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_Cerebellum_FT_70dpf_4_CAAGCCAA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_70dpf_4.bed

## Cerebellum NB

sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_M_1_GTGATCCA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_1.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_M_2_TCCAACTG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_2.bed
sed '1,11d' ../RRBS/SSR_INRA_GS_WP1_cerebellum_NB_F_3_AAGGACCA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_cerebellum_NB_3.bed

## Skin 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_1_AACCGAAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_3_GTACCACA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Skin_FT_30dpf_POOL_6_AACCGAAC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_30dpf_POOL_6.bed

## Skin 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_1_TCACTCGA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_2_TCCTGGTA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Back_skin_FT_70dpf_3_AGGTTCCT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_70dpf_3.bed

## Skin NB

sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_F_3_CTATCCAC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_3.bed
sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_M_1_CTATCCAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_back_skin_NB_M_2_ACTGCACT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_2.bed

## Lung 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_2_CTCGAACA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_3_GTCTCATC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_30dpf_POOL_6_CTCGAACA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_30dpf_POOL_6.bed

## Lung 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_1_TCCATTGC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_2_AACACTGG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_FT_70dpf_3_TCCATTGC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_lung_70dpf_3.bed

## Lung NB

sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_F_4_ATAACGCC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_M_1_ATAACGCC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Lung_NB_M_2_TCTAGTCC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_skin_NB_2.bed

## Muscle 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_2_CCTTAGGT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_3_ACGATCAG_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Hindlimb_muscle_FT_30dpf_POOL_6_CCTTAGGT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_30dpf_POOL_6.bed

## Muscle 70dpf

sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_1_AAGTCCTC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_2_AGCAGACA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_hindlimb_muscle_FT_70dpf_4_AAGTCCTC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_70dpf_4.bed

## Muscle NB

sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_F_4_ACGCTTCT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_M_1_ACGCTTCT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_gluteus_medius_NB_M_2_GACATCTC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_muscle_NB_2.bed

## Kidney 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_1_ACGGACTT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_3_GCCAGAAT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_30dpf_POOL_6_ACGGACTT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_30dpf_POOL_6.bed

## Kidney 70dpf

sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_1_CTGACTAC_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_2_CATCAACC_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_Kidney_FT_70dpf_3_GAACCTTC_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_70dpf_3.bed

## Kidney NB

sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_F_3_CCGGAATA_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_3.bed
sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_M_1_CCGGAATA_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_kidney_NB_M_2_CTCGACTT_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_kidney_NB_2.bed

## Intestine 30dpf

sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_1_TCTAGGAG_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_3_ACTCTCCA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_3.bed
sed '1,11d' SSR_INRA_GS_WP1_Small_Intestine_FT_30dpf_POOL_6_TCTAGGAG_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_30dpf_POOL_6.bed

## Intestine 70dpf

sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_1_AACGCCTT_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_1.bed
sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_2_CGCAACTA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_2.bed
sed '1,11d' SSR_INRA_GS_WP1_small_intestine_FT_70dpf_4_AACGCCTT_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_70dpf_4.bed

## Intestine NB

sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_F_4_CCAAGTAG_L003.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_4.bed
sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_M_1_CCAAGTAG_L001.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_1.bed
sed '1,11d' SSR_INRA_GS_WP1_Ileum_NB_M_2_CTAGCTCA_L002.f.bed | sed 's/chr//g' | awk '{OFS="\t"} {if ($5=="CG") print $0}' > CpG_deg_intestine_NB_2.bed
((current_step++))

# Step 6: Finalize formatting

progress "Step 6: Finalize formatting" "$current_step" "$total_steps"
echo -e "chromosome\tgene\tpval_1-2\tadjpval_1-2\tlogFC_1-2\tlogCPM_1-2\tpval_2-3\tadjpval_2-3\tlogFC_2-3\tlogCPM_2-3\tmeth_1\tmeth_2\tmeth_3" > header.txt

## Liver 30dpf

bedtools intersect -wo -a CpG_deg_liver_30dpf_POOL_2.bed -b CpG_deg_liver_30dpf_POOL_3.bed > ISEC_liver_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_liver_30dpf_2_3_temp.bed -b  CpG_deg_liver_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_liver_30dpf.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_POOLS_merged_liver_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_liver_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_liver_30dpf.bed > final_liver_30dpf.bed

## Liver 70dpf

bedtools intersect -wo -a CpG_deg_liver_70dpf_1.bed -b CpG_deg_liver_70dpf_2.bed > ISEC_liver_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_liver_70dpf_1_2_temp.bed -b  CpG_deg_liver_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_liver_70dpf.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_merged_liver_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_liver_70dpf.bed
cat header.txt temp_CpG_genes_merged_liver_70dpf.bed > final_liver_70dpf.bed

## Liver NB

bedtools intersect -wo -a CpG_deg_liver_NB_1.bed -b CpG_deg_liver_NB_2.bed > ISEC_liver_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_liver_NB_1_2_temp.bed -b  CpG_deg_liver_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_liver_NB.bed
bedtools intersect -wo -a ISEC_deg_liver_1_2_3.bed -b CpG_merged_liver_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_liver_NB.bed
cat header.txt temp_CpG_genes_merged_liver_NB.bed > final_liver_NB.bed

## Hindbrain 30dpf

bedtools intersect -wo -a CpG_deg_hindbrain_30dpf_POOL_2.bed -b CpG_deg_hindbrain_30dpf_POOL_3.bed > ISEC_hindbrain_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_hindbrain_30dpf_2_3_temp.bed -b  CpG_deg_hindbrain_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_hindbrain_30dpf.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_POOLS_merged_hindbrain_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_hindbrain_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_hindbrain_30dpf.bed > final_hindbrain_30dpf.bed

## Hindbrain 70dpf

bedtools intersect -wo -a CpG_deg_cerebellum_70dpf_1.bed -b CpG_deg_cerebellum_70dpf_2.bed > ISEC_cerebellum_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_cerebellum_70dpf_1_2_temp.bed -b  CpG_deg_cerebellum_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_cerebellum_70dpf.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_merged_cerebellum_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_cerebellum_70dpf.bed
cat header.txt temp_CpG_genes_merged_cerebellum_70dpf.bed > final_hindbrain_70dpf.bed

## Hindbrain NB

bedtools intersect -wo -a CpG_deg_cerebellum_NB_1.bed -b CpG_deg_cerebellum_NB_2.bed > ISEC_cerebellum_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_cerebellum_NB_1_2_temp.bed -b  CpG_deg_cerebellum_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_cerebellum_NB.bed
bedtools intersect -wo -a ISEC_deg_cerebellum_1_2_3.bed -b CpG_merged_cerebellum_NB.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_cerebellum_NB.bed
cat header.txt temp_CpG_genes_merged_cerebellum_NB.bed > final_hindbrain_NB.bed

## Muscle 30dpf

bedtools intersect -wo -a CpG_deg_muscle_30dpf_POOL_2.bed -b CpG_deg_muscle_30dpf_POOL_3.bed > ISEC_muscle_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_muscle_30dpf_2_3_temp.bed -b  CpG_deg_muscle_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_muscle_30dpf.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_POOLS_merged_muscle_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_muscle_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_muscle_30dpf.bed > final_muscle_30dpf.bed

## Muscle 70dpf

bedtools intersect -wo -a CpG_deg_muscle_70dpf_1.bed -b CpG_deg_muscle_70dpf_2.bed > ISEC_muscle_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_muscle_70dpf_1_2_temp.bed -b  CpG_deg_muscle_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_muscle_70dpf.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_merged_muscle_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_muscle_70dpf.bed
cat header.txt temp_CpG_genes_merged_muscle_70dpf.bed > final_muscle_70dpf.bed

## Muscle NB

bedtools intersect -wo -a CpG_deg_muscle_NB_1.bed -b CpG_deg_muscle_NB_2.bed > ISEC_muscle_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_muscle_NB_1_2_temp.bed -b  CpG_deg_muscle_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_muscle_NB.bed
bedtools intersect -wo -a ISEC_deg_muscle_1_2_3.bed -b CpG_merged_muscle_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_muscle_NB.bed
cat header.txt temp_CpG_genes_merged_muscle_NB.bed > final_muscle_NB.bed

## Kidney 30dpf

bedtools intersect -wo -a CpG_deg_kidney_30dpf_POOL_1.bed -b CpG_deg_kidney_30dpf_POOL_3.bed > ISEC_kidney_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_kidney_30dpf_1_3_temp.bed -b  CpG_deg_kidney_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_kidney_30dpf.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_POOLS_merged_kidney_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_kidney_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_kidney_30dpf.bed > final_kidney_30dpf.bed

## Kidney 70dpf

bedtools intersect -wo -a CpG_deg_kidney_70dpf_1.bed -b CpG_deg_kidney_70dpf_2.bed > ISEC_kidney_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_kidney_70dpf_1_2_temp.bed -b  CpG_deg_kidney_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_kidney_70dpf.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_merged_kidney_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_kidney_70dpf.bed
cat header.txt temp_CpG_genes_merged_kidney_70dpf.bed > final_kidney_70dpf.bed

## Kidney NB

bedtools intersect -wo -a CpG_deg_kidney_NB_1.bed -b CpG_deg_kidney_NB_2.bed > ISEC_kidney_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_kidney_NB_1_2_temp.bed -b  CpG_deg_kidney_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_kidney_NB.bed
bedtools intersect -wo -a ISEC_deg_kidney_1_2_3.bed -b CpG_merged_kidney_NB.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_kidney_NB.bed
cat header.txt temp_CpG_genes_merged_kidney_NB.bed > final_kidney_NB.bed

## Lung 30dpf

bedtools intersect -wo -a CpG_deg_lung_30dpf_POOL_2.bed -b CpG_deg_lung_30dpf_POOL_3.bed > ISEC_lung_30dpf_2_3_temp.bed
bedtools intersect -wo -a ISEC_lung_30dpf_2_3_temp.bed -b  CpG_deg_lung_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_lung_30dpf.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_POOLS_merged_lung_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_lung_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_lung_30dpf.bed > final_lung_30dpf.bed

## Lung 70dpf

bedtools intersect -wo -a CpG_deg_lung_70dpf_1.bed -b CpG_deg_lung_70dpf_2.bed > ISEC_lung_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_lung_70dpf_1_2_temp.bed -b  CpG_deg_lung_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_lung_70dpf.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_merged_lung_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_lung_70dpf.bed
cat header.txt temp_CpG_genes_merged_lung_70dpf.bed > final_lung_70dpf.bed

## Lung NB

bedtools intersect -wo -a CpG_deg_skin_NB_1.bed -b CpG_deg_skin_NB_2.bed > ISEC_lung_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_lung_NB_1_2_temp.bed -b  CpG_deg_skin_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_lung_NB.bed
bedtools intersect -wo -a ISEC_deg_lung_1_2_3.bed -b CpG_merged_lung_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_lung_NB.bed
cat header.txt temp_CpG_genes_merged_lung_NB.bed > final_lung_NB.bed

## Skin 30dpf

bedtools intersect -wo -a CpG_deg_skin_30dpf_POOL_1.bed -b CpG_deg_skin_30dpf_POOL_3.bed > ISEC_skin_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_skin_30dpf_1_3_temp.bed -b  CpG_deg_skin_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_skin_30dpf.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_POOLS_merged_skin_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_skin_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_skin_30dpf.bed > final_skin_30dpf.bed

## Skin 70dpf

bedtools intersect -wo -a CpG_deg_skin_70dpf_1.bed -b CpG_deg_skin_70dpf_2.bed > ISEC_skin_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_skin_70dpf_1_2_temp.bed -b  CpG_deg_skin_70dpf_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_skin_70dpf.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_merged_skin_70dpf.bed| awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_skin_70dpf.bed
cat header.txt temp_CpG_genes_merged_skin_70dpf.bed > final_skin_70dpf.bed

## Skin NB

bedtools intersect -wo -a CpG_deg_skin_NB_1.bed -b CpG_deg_skin_NB_2.bed > ISEC_skin_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_skin_NB_1_2_temp.bed -b  CpG_deg_skin_NB_3.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_skin_NB.bed
bedtools intersect -wo -a ISEC_deg_skin_1_2_3.bed -b CpG_merged_skin_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_skin_NB.bed
cat header.txt temp_CpG_genes_merged_skin_NB.bed > final_skin_NB.bed

## Intestine 30dpf

bedtools intersect -wo -a CpG_deg_intestine_30dpf_POOL_1.bed -b CpG_deg_intestine_30dpf_POOL_3.bed > ISEC_intestine_30dpf_1_3_temp.bed
bedtools intersect -wo -a ISEC_intestine_30dpf_1_3_temp.bed -b  CpG_deg_intestine_30dpf_POOL_6.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_POOLS_merged_intestine_30dpf.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_POOLS_merged_intestine_30dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_POOLs_genes_merged_intestine_30dpf.bed
cat header.txt temp_CpG_POOLs_genes_merged_intestine_30dpf.bed > final_intestine_30dpf.bed

## Intestine 70dpf

bedtools intersect -wo -a CpG_deg_intestine_70dpf_1.bed -b CpG_deg_intestine_70dpf_2.bed > ISEC_intestine_70dpf_1_2_temp.bed
bedtools intersect -wo -a ISEC_intestine_70dpf_1_2_temp.bed -b  CpG_deg_intestine_70dpf_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_intestine_70dpf.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_merged_intestine_70dpf.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_intestine_70dpf.bed
cat header.txt temp_CpG_genes_merged_intestine_70dpf.bed > final_intestine_70dpf.bed

## Intestine NB

bedtools intersect -wo -a CpG_deg_intestine_NB_1.bed -b CpG_deg_intestine_NB_2.bed > ISEC_intestine_NB_1_2_temp.bed
bedtools intersect -wo -a ISEC_intestine_NB_1_2_temp.bed -b  CpG_deg_intestine_NB_4.bed | awk '{OFS="\t"} {print $1,$2,$3,$7,$16,$26}' > CpG_merged_intestine_NB.bed
bedtools intersect -wo -a ISEC_deg_ileum_1_2_3.bed -b CpG_merged_intestine_NB.bed | awk '{OFS="\t"} {print $1,$4,$5,$6,$7,$8,$9,$10,$11,$12,$16,$17,$18}' > temp_CpG_genes_merged_intestine_NB.bed
cat header.txt temp_CpG_genes_merged_intestine_NB.bed > final_intestine_NB.bed

# Finally, update the progress bar to 100% when the script is done.
progress "Progress: Script completed" "$total_steps" "$total_steps"
exit 0