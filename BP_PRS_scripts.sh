
# Target data quality control
## Most of them have been described in https://github.com/neurogenetics/GWAS-pipeline
## Soft call imputed genotyping data was used in PRS calculation
eval mkdir soft_calls
for chr in {1..22};do # do this for chr1 to chr22
    gunzip raw/chr$chr.info.gz # files from imputation server
    plink --vcf raw/chr$chr.dose.vcf.gz --make-bed --out s1_chr$chr --const-fid
    plink --bfile s1_chr$chr --bmerge s1_chr$chr --merge-mode 6 --out s8_chr$chr
    plink --bfile s1_chr$chr --exclude s8_chr$chr.diff --make-bed --out s2_chr$chr
    plink --bfile s2_chr$chr --list-duplicate-vars --out s4_chr$chr
    plink --bfile s2_chr$chr --exclude s4_chr$chr.dupvar --make-bed --out softcalls_chr$chr
    plink --bfile softcalls_chr$chr --qual-scores raw/chr$chr.info 7 1 1 --qual-threshold 0.3 --make-bed --out soft_calls/${NAME}_chr$chr
done
rm s?_chr$chr*
rm softcalls_chr$chr*

# merge all chromosomes using PLINK
ls soft_calls/*.bed | sed 's/\.bed//' > soft_calls/merge_list
plink --merge-list soft_calls/merge_list --make-bed --out soft_calls/ABCD.all.soft

# Do basic QC on the merged data
# reference: https://choishingwan.github.io/PRS-Tutorial/
plink \
    --bfile soft_calls/ABCD.all.soft \
    --maf 0.01 \
    --hwe 1e-6 \
    --geno 0.01 \
    --out soft_calls/ABCD.all.soft.qc
# exclude duplicated snps using PLINK2
plink2 \
    --bfile soft_calls/ABCD.all.soft.qc \
    --rm-dup force-first \
    --make-bed \
    --out soft_calls/ABCD.all.soft.qc.nodup

# PRS calculation for ABCD dataset using PGC bipolar summary statistics
# use 0.1 as the pvalue threshold which is specified in the paper.
base=pgc-bip2021-all.vcf.tsv
target=soft_calls/ABCD.all.soft.qc.nodup
output_folder=/path/to/your/output/folder
output_name=ABCD_BP_PRSice
out=${output_folder}/${output_name}
Rscript /PRSice/PRSice.R \
    --prsice /PRSice/PRSice_linux \
    --a1 A1 \
    --a2 A2 \
    --no-full \
    --bar-levels 0.1 \
    --base $base \
    --beta  \
    --fastscore \
    --no-regress \
    --print-snp  \
    --clump-kb 250kb \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --num-auto 22 \
    --out $out \
    --pvalue PVAL \
    --score avg \
    --snp ID \
    --stat BETA \
    --target $target \
    --ultra 
