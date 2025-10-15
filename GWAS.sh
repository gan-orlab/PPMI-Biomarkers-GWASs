### GWAS script for PPMI biomarkers
## Includes both original PD/HC/SWEDD cohort (pdhc) and prodromal/genetic cohort biomarkers (progen/prgn)

## ==== Set up directories ====
# run this in your working directory
mkdir pdhc 
mkdir pdhc/plots
mkdir pdhc/sumstats
mkdir pdhc/linear
mkdir pdhc/annos
mkdir pdhc/distribution_data
mkdir progen
mkdir progen/plots
mkdir progen/sumstats
mkdir progen/linear
mkdir progen/annos
mkdir progen/distribution_data
mkdir meta
mkdir meta/linear
mkdir meta/sumstats 
mkdir meta/plots 
mkdir meta/annos

# replace below with wherever your info files are
cp ~/runs/go_lab/GRCh37/PPMI/pdhc_pdhc/imputed_topmed/pdhc.info pdhc.info
cp ~/runs/go_lab/GRCh37/PPMI/prodromal_genetic/imputed_topmed/prgn.info prgn.info

## ==== GWAS for biomarkers in the PD/HC/SWEDD dataset ====
# phenos_pdhc.txt is a file listing the phenotypes of biomarkers with data in the pdhc cohort
cat phenos/phenos_pdhc.txt | while read line
do 
    ## Distribution data and plots
    echo ${line} > marker.txt
    awk '{print $1}' ppmi_hg38_pdhc_noindels.fam > pdhc_grep.txt
    awk '{print $1}' prgn_hg38_sc_pdhc_noindels.fam > progen_grep.txt
    grep -f pdhc_grep.txt phenos_new/${line}.txt > biomarker_data.txt
    R < distribution_temp.R --no-save
    mv bio_hist.tiff pdhc/distribution_data/hist_pdhc_${line}.tiff
    paste shapiro.txt >> pdhc/distribution_data/shapiro_pdhc_all.txt
 
    ## Linear regression 
    plink --bfile ppmi_hg38_pdhc_noindels --pheno phenos/${line}.txt --update-sex pdhc_sex.txt --make-bed --out pdhc_${line}
    plink --bfile pdhc_${line} --linear --ci .95 --maf 0.05 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out pdhc_sc_tm_maf05_${line}
    awk '{if ($12!="NA") print $0}' pdhc_sc_tm_maf05_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.pdhc_sc_tm_maf05_${line}.assoc.linear

    # Additional line to manually add GBA1 N370S (MAF would normally be below threshold)
    plink --bfile ppmi_hg38_pdhc_noindels --pheno phenos/${line}.txt --update-sex  pdhc_sex.txt --extract GBA_N370S.txt --make-bed --out pdhc_GBA_${line}
    plink --bfile pdhc_GBA_${line} --linear --ci .95 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out pdhc_sc_tm_maf05_GBA_${line}
    awk '{if ($12!="NA") print $0}' pdhc_sc_tm_maf05_GBA_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.pdhc_sc_tm_maf05_GBA_${line}.assoc.linear
    cat p.pdhc_sc_tm_maf05_GBA_${line}.assoc.linear p.pdhc_sc_tm_maf05_${line}.assoc.linear > pdhc/linear/p.pdhc_sc_tm_maf05_${line}.assoc.linear
    cat pdhc_sc_tm_maf05_GBA_${line}.assoc.linear pdhc_sc_tm_maf05_${line}.assoc.linear > pdhc/linear/pdhc_sc_tm_maf05_${line}.assoc.linear

    # Clean up 
    mv *.log log/
    rm pdhc_${line}.*
    rm pdhc_GBA_${line}.*
    rm p.* 
    rm pdhc_sc_tm_maf05_*

    # Generate summary statistics and plots
    echo ${line} > marker.txt 
    cp pdhc.info info
    cp pdhc/linear/p.pdhc_sc_tm_maf05_${line}.assoc.linear assoc
    R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save #find this script at https://github.com/gan-orlab/GlcCer_GWAS/blob/main/sumstats_from_assoc.R
    paste lambda.txt >> pdhc/plots/pdhc_lambdas.txt
    mv QQ.tiff pdhc/plots/QQ_pdhc_sc_tm_maf05_${line}.tiff
    mv ManH.tiff pdhc/plots/ManH_pdhc_sc_tm_maf05_${line}.tiff
    mv Metal.tab pdhc/sumstats/Metal_pdhc_sc_tm_maf05_${line}.txt
    mv COJO.tab pdhc/sumstats/COJO_pdhc_sc_tm_maf05_${line}.txt
    mv fullSTATS.tab pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt

    # Pull annotations for summary statistics using ANNOVAR
    head -n1 pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt > header.txt
    awk '{print $2,$2,$2,$5,$4}' pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt | sed '1d' > input.annovar.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
        -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210502,gnomad_genome \
        -operation g,f,f,f -nastring . -polish
    paste pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt annotatedtable.hg38_multianno.txt > temp.txt
    tr ' ' '\t' < temp.txt > temp2.txt
    sort -gk10 temp2.txt > pdhc/annos/ANNO_pdhc_sc_tm_maf05_${line}.txt
done

## ==== GWAS for biomarkers in both PD/HC/SWEDD and prodromal/genetic datasets ====
# phenos_meta.txt is a file listing the phenotypes of biomarkers with data in both cohorts
cat phenos/phenos_meta.txt | while read line
do 
    ## Distribution data and plots
    # Prep work
    echo ${line} > marker.txt
    awk '{print $1}' ppmi_hg38_pdhc_noindels.fam > pdhc_grep.txt
    awk '{print $1}' prgn_hg38_sc_pdhc_noindels.fam > progen_grep.txt

    # Run for pdhc
    grep -f pdhc_grep.txt phenos_new/${line}.txt > biomarker_data.txt
    R < distribution_temp.R --no-save
    mv bio_hist.tiff pdhc/distribution_data/hist_pdhc_${line}.tiff
    paste shapiro.txt >> pdhc/distribution_data/shapiro_pdhc_all.txt

    # Run for progen
    grep -f progen_grep.txt phenos_new/${line}.txt > biomarker_data.txt
    R < distribution_temp.R --no-save
    mv bio_hist.tiff progen/distribution_data/hist_prgn_${line}.tiff
    paste shapiro.txt >> progen/distribution_data/shapiro_prgn_all.txt

    ## Linear regression 
    # pdhc
    plink --bfile ppmi_hg38_pdhc_noindels --pheno phenos/${line}.txt --update-sex pdhc_sex.txt --make-bed --out pdhc_${line}
    plink --bfile pdhc_${line} --linear --ci .95 --maf 0.05 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out pdhc_sc_tm_maf05_${line}
    awk '{if ($12!="NA") print $0}' pdhc_sc_tm_maf05_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.pdhc_sc_tm_maf05_${line}.assoc.linear
    plink --bfile ppmi_hg38_pdhc_noindels --pheno phenos_new/${line}.txt --update-sex pdhc_sex.txt --extract GBA_N370S.txt --make-bed --out pdhc_GBA_${line}
    plink --bfile pdhc_GBA_${line} --linear --ci .95 --covar pdhc_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out pdhc_sc_tm_maf05_GBA_${line}
    awk '{if ($12!="NA") print $0}' pdhc_sc_tm_maf05_GBA_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.pdhc_sc_tm_maf05_GBA_${line}.assoc.linear
    cat p.pdhc_sc_tm_maf05_GBA_${line}.assoc.linear p.pdhc_sc_tm_maf05_${line}.assoc.linear > pdhc/linear/p.pdhc_sc_tm_maf05_${line}.assoc.linear
    cat pdhc_sc_tm_maf05_GBA_${line}.assoc.linear pdhc_sc_tm_maf05_${line}.assoc.linear > pdhc/linear/pdhc_sc_tm_maf05_${line}.assoc.linear

    # prgn
    plink --bfile prgn_hg38_sc_pdhc_noindels --pheno phenos/${line}.txt --update-sex prgn_sex.txt --make-bed --out prgn_${line}
    plink --bfile prgn_${line} --linear --ci .95 --maf 0.05 --covar prgn_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out prgn_sc_tm_maf05_${line}
    awk '{if ($12!="NA") print $0}' prgn_sc_tm_maf05_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.prgn_sc_tm_maf05_${line}.assoc.linear
    plink --bfile prgn_hg38_sc_pdhc_noindels --pheno phenos/${line}.txt --update-sex prgn_sex.txt --extract GBA_N370S.txt --make-bed --out prgn_GBA_${line}
    plink --bfile prgn_GBA_${line} --linear --ci .95 --maf 0.05 --covar prgn_covar_2.txt --covar-name Sex,Age,Status,PC1,PC2,PC3,PC4,PC5 --out prgn_sc_tm_maf05_GBA_${line}
    awk '{if ($12!="NA") print $0}' prgn_sc_tm_maf05_GBA_${line}.assoc.linear | grep 'ADD' | sort -gk12 > p.prgn_sc_tm_maf05_GBA_${line}.assoc.linear
    cat p.prgn_sc_tm_maf05_GBA_${line}.assoc.linear p.prgn_sc_tm_maf05_${line}.assoc.linear > progen/linear/p.prgn_sc_tm_maf05_${line}.assoc.linear
    cat prgn_sc_tm_maf05_GBA_${line}.assoc.linear prgn_sc_tm_maf05_${line}.assoc.linear > progen/linear/prgn_sc_tm_maf05_${line}.assoc.linear

    # Clean up 
    rm pdhc_${line}.*
    rm pdhc_GBA_${line}.*
    rm p.* 
    rm pdhc_sc_tm_maf05_*
    rm prgn_${line}.*
    rm prgn_GBA_${line}.*
    rm p.*
    rm prgn_sc_tm_maf05*

    # Generate summary statistics and plots
    # pdhc
    echo ${line} > marker.txt 
    cp pdhc.info info
    cp pdhc/linear/p.pdhc_sc_tm_maf05_${line}.assoc.linear assoc
    R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
    paste lambda.txt >> pdhc/plots/pdhc_lambdas.txt
    mv QQ.tiff pdhc/plots/QQ_pdhc_sc_tm_maf05_${line}.tiff
    mv ManH.tiff pdhc/plots/ManH_pdhc_sc_tm_maf05_${line}.tiff
    mv Metal.tab pdhc/sumstats/Metal_pdhc_sc_tm_maf05_${line}.txt
    mv COJO.tab pdhc/sumstats/COJO_pdhc_sc_tm_maf05_${line}.txt
    mv fullSTATS.tab pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt

    # prgn
    echo ${line} > marker.txt     
    cp prgn.info info
    cp progen/linear/p.prgn_sc_tm_maf05_${line}.assoc.linear assoc
    R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
    paste lambda.txt >> progen/plots/prgn_lambdas.txt
    mv QQ.tiff progen/plots/QQ_prgn_sc_tm_maf05_${line}.tiff
    mv ManH.tiff progen/plots/ManH_prgn_sc_tm_maf05_${line}.tiff
    mv Metal.tab progen/sumstats/Metal_prgn_sc_tm_maf05_${line}.txt
    mv COJO.tab progen/sumstats/COJO_prgn_sc_tm_maf05_${line}.txt
    mv fullSTATS.tab progen/sumstats/fullSTATS_prgn_sc_tm_maf05_${line}.txt

    # Clean up 
    rm assoc 
    rm info 
    rm marker.txt 
    rm lambda.txt

    # Pull annotations for summary statistics using ANNOVAR
    # pdhc
    head -n1 pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt > header.txt
    awk '{print $1,$2,$2,$5,$4}' pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt | sed '1d' > input.annovar.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
        -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210502,gnomad_genome  \
        -operation g,f,f,f -nastring . -polish
    paste pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt annotatedtable.hg38_multianno.txt > temp.txt
    tr ' ' '\t' < temp.txt > temp2.txt
    sort -gk10 temp2.txt > pdhc/annos/ANNO_pdhc_sc_tm_maf05_${line}.txt

    # prgn
    head -n1 progen/sumstats/fullSTATS_prgn_sc_tm_maf05_${line}.txt > header.txt
    awk '{print $1,$2,$2,$5,$4}' progen/sumstats/fullSTATS_prgn_sc_tm_maf05_${line}.txt | sed '1d' > input.annovar.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
        -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210502,gnomad_genome  \
        -operation g,f,f,f -nastring . -polish
    paste progen/sumstats/fullSTATS_prgn_sc_tm_maf05_${line}.txt annotatedtable.hg38_multianno.txt > temp.txt
    tr ' ' '\t' < temp.txt > temp2.txt
    sort -gk10 temp2.txt > progen/annos/ANNO_prgn_sc_tm_maf05_${line}.txt

    # Clean up
    rm header.txt 
    rm annotatedtable.hg38_multianno.txt 
    rm input.annovar.txt
done 

## ==== Meta-analysis for biomarkers in both datasets ====
cat phenos/phenos_meta.txt | while read line
do 
    ## Get sumstats in workable variables for R script
    cp pdhc/sumstats/fullSTATS_pdhc_sc_tm_maf05_${line}.txt meta1.tab
    cp progen/sumstats/fullSTATS_prgn_sc_tm_maf05_${line}.txt meta2.tab 

    # Flip mismatch alleles in pdhc summary stats
    R < allele_flip.R --no-save
    awk '{print $1,$2,$1":"$2":"$4":"$5,$5,$4,$14,$7,$8,$15,$10,$11,$12,$13}' mismatches.txt > temp.txt
    tr ' ' '\t' < temp.txt > temp2.txt
    sed '1d' temp2.txt > temp.txt 
    cat matches.txt temp.txt > meta1.tab

    # Run metal and format output
    module load StdEnv/2020
    module load gcc/9.3.0
    module load metal
    metal ~/runs/emsom/scripts/metal.txt 
    awk 'BEGIN{FS=OFS="\t"}{split($1,snp,":"); print snp[1],snp[2]}' metal1.tbl > chrbp.txt
    sed 's/chr//g' chrbp.txt > chrbp2.txt
    paste chrbp2.txt metal1.tbl > meta/linear/meta_ppmi_sc_tm_maf05_${line}.tbl 

    # Clean up
    rm chrbp.txt
    rm metal1.tbl
    rm meta1.tab
    rm meta2.tab
    rm temp.*

    # Generate summary statistics and plots
    cp meta/linear/meta_ppmi_sc_tm_maf05_${line}.tbl meta.tab
    R < ~/runs/emsom/scripts/sumstats_from_meta.R --no-save
    mv fullSTATS.meta.tab meta/sumstats/fullSTATS_meta_ppmi_sc_tm_maf05_${line}.txt
    mv QQ.meta.tiff meta/plots/QQ_meta_ppmi_sc_tm_maf05_${line}.tiff
    mv ManH.meta.tiff meta/plots/ManH_meta_ppmi_sc_tm_maf05_${line}.tiff

    # Pull annotations for summary statistics using ANNOVAR
    head -n1 meta/sumstats/fullSTATS_meta_ppmi_sc_tm_maf05_${line}.txt > header.txt 
    awk '{print $1,$2,$2,$5,$4}' meta/sumstats/fullSTATS_meta_ppmi_sc_tm_maf05_${line}.txt | sed '1d' > input.annovar.txt
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
        -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210502,gnomad_genome  \
        -operation g,f,f,f -nastring . -polish
    paste meta/sumstats/fullSTATS_meta_ppmi_sc_tm_maf05_${line}.txt annotatedtable.hg38_multianno.txt > temp.txt
    tr ' ' '\t' < temp.txt > temp2.txt
    sort -gk10 temp2.txt > meta/annos/ANNO_meta_ppmi_sc_tm_maf05_${line}.txt

    # Clean up
    rm header.txt 
    rm annotatedtable.hg38_multianno.txt 
    rm input.annovar.txt
done

