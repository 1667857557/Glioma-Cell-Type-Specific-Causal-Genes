###s-LDSC------
source activate ldsc
export PATH=$PATH:/home/huang/anaconda3/bin
cd /mnt/d/ldsc

for file in /mnt/d/ldsc/corces/*.annot.gz; do
input=$(basename "$file" | cut -d. -f1-$(($(basename "$file" | grep -o "\." | wc -l) - 1)))
digit=$(echo "$input" | grep -oE "[0-9]+$")
./ldsc.py \
--l2 \
--bfile /mnt/d/ldsc/1000G_EUR_Phase3_plink/1000G.EUR.QC.$digit \
--ld-wind-cm 1 \
--annot "$file" \
--print-snps /mnt/d/ldsc/LDSCORE_list.txt \
--thin-annot \
--out "/mnt/d/ldsc/corces/$input"
done

for CHR in {1..22}
do
./make_annot.py \
--bed-file /mnt/d/ldsc/corces.bed \
--bimfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.${CHR}.bim \
--annot-file corces.${CHR}.annot.gz
done

directory_path="/mnt/g/Linux/ldsc/Kim_ATAC/annots/corces"
seen_names=()
for full_filepath in "$directory_path"/*.annot.gz; do
full_filename=$(basename "$full_filepath")
object_name=$(echo "$full_filename" | cut -d'.' -f1)
if [ -f "$full_filepath" ]; then
if [[ " ${seen_names[@]} " =~ " $object_name " ]]; then
continue
fi
seen_names+=("$object_name")

./ldsc.py \
--h2 GBM_TCSC.sumstats \
--ref-ld-chr /mnt/g/Linux/ldsc/Kim_ATAC/annots/corces/${object_name}.,/mnt/g/Linux/ldsc/baselineLD_v2.3/baselineLD. \
--w-ld-chr /mnt/g/Linux/ldsc/weights_hm3_no_hla/weights. \
--overlap-annot \
--frqfile-chr /mnt/g/Linux/ldsc/frq/1000G.EUR.QC. \
--print-coefficients \
--print-delete-vals \
--out "/mnt/g/Linux/ldsc/Kim_ATAC/annots/corces/${object_name}."
fi
done
