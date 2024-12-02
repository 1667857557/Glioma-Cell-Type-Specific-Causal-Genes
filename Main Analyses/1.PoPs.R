###Whole tissue features pops----
export PATH=$PATH:/mnt/data/userdata/tyu-bf0hxdy3yggq/pops1

###MAGMA
./magma \
--bfile ./data/1000G.EUR \
--gene-annot ./data/magma_0kb.genes.annot \
--pval nonGBM.txt ncol=N \
--gene-model snp-wise=mean \
--out nonGBM

###PoPs

python3.10 ./pops.feature_selection.py \
--features ./data/PoPS.features.txt.gz \
--gene_results ./data/ALL_GLIOMA \
--out ALL_GLIOMA

for CHR in {1..22}
do
python3.10 ./pops.predict_scores.py \
--gene_loc ./data/gene_loc.txt \
--gene_results ./data/ALL_GLIOMA \
--features ./data/PoPS.features.txt.gz \
--selected_features ALL_GLIOMA.features \
--control_features ./data/control.features \
--chromosome ${CHR} \
--out ALL_GLIOMA${CHR}
done

####Brain-specific features POPS----
setwd("D:/linux/pops1/data/features_munged/")
B <- vroom("D:/linux/pops1/data/features_munged/all/full_features_jul17.txt")
colnames(gwas.sumstats)
A<-vroom("brain_gene_features_metadata.txt")
C <- B[, c("ENSGID",intersect(names(B), A$...1))]
write_tsv(C,file =  "brain_features_jul17.txt")

####POPS
export PATH=$PATH:/mnt/g/linux/FLAMES/

python3 pops.py \
--gene_annot_path /mnt/g/linux/FLAMES/pops_features_full/gene_annot.txt \
--feature_mat_prefix /mnt/g/linux/FLAMES/pops_features_full/munged_features/pops_features \
--num_feature_chunks 99 \
--magma_prefix /mnt/g/linux/FLAMES/GBM \
--control_features /mnt/g/linux/FLAMES/pops_features_full/control.features \
--out_prefix /mnt/g/linux/FLAMES/GBM

python3 munge_feature_directory.py \
--gene_annot_path /mnt/g/linux/FLAMES/gene_annot_jun10.txt \
--feature_dir /mnt/g/linux/FLAMES/POPS_data/features_brain_munged/ \
--save_prefix /mnt/g/linux/FLAMES/POPS_data/features_brain_munged/pops_brain_features \
--max_cols 5000