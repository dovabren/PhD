bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE[\t%GT]\n' calc_and_other_ceratina_filter3.vcf > calc_and_other_ceratina_filter3_GT.vcf

        #remove 14 and 16 

        awk '{$14=$16=""; print $0}' calc_and_other_ceratina_filter3_GT.vcf > calc_and_other_ceratina_filter3_GT_rm_highalt.vcf
        
#try keeping last 2 samples
samples.txt -- add in Sample_19 Sample 20


sed -i 's/ /\t/g' calc_and_other_ceratina_filter3_GT_rm_highalt.vcf
sed -i 's/\t\t/\t/g' calc_and_other_ceratina_filter3_GT_rm_highalt.vcf

awk 'BEGIN {IFS = OFS = "\t"} {print $1,$2,".",$3,$4,".",$5,".","GT",$0}' calc_and_other_ceratina_filter3_GT_rm_highalt.vcf| cut -f 1-9,15- | sed '1i#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample_1\tSample_2\tSample_3\tSample_4\tSample_5\tSample_6\tSample_7\tSample_8\tSample_9\tSample_10\tSample_11\tSample_12\tSample_13\tSample_14\tSample_15\tSample_16\tSample_17\tSample_18\tOut_1\tOut_2' | sed '1i##fileformat=VCFv4.2' > Ceratina.vcf

#Keep samples with 80% coverage and max two alleles ## july 28 2022- change max missing to 0.7
vcftools --vcf Ceratina.vcf --max-missing 0.7 --keep samples.txt --max-alleles 2 --recode --out TargetFilter
#After filtering, kept 48217 out of a possible 59628 Sites

#Outgroup filtering if fixed, only one needed, if polymorphic both needed.
awk 'BEGIN {IFS = OFS = "\t"} {if ($30 == "0/0" || $30 == "1/1"|| $31 == "0/0" || $31 == "1/1") print $1,$2,$30,$31; else if ($30 == "0/1" && $31 == "0/1") print $1,$2,$30,$31; else print $1,$2,$30,$31,"Check"}' Ceratina.vcf | sed '/Check/d' > OutgroupFilter.txt

#Merge
sed '/#/d' TargetFilter.recode.vcf | cut -f 1,2 | sed 's/\t/_/' > TargetPos.txt
cut -f 1,2 OutgroupFilter.txt | sed 's/\t/_/' > OutPos.txt

R
setwd('/data3/african/sideproject/dova')

tar = read.table(file='TargetPos.txt', header=F)
out = read.table(file='OutPos.txt', header=F)
m = merge(tar, out, by="V1")
write.table(m, "Matches.txt", sep='\t', quote=F, row.names=F)

sed -i 's/_/\t/2; 1d' Matches.txt 


#email katie about files jul 28
#Select longest gene
sort -nrk 3 -k1 mrnaselection.txt > sorted
awk '! a[$1]++' sorted > Genes.txt


#Extract samples within exons
bedtools intersect -a intervals.txt -b Ceratina.vcf -wa -wb | cut -f 4,5,6 > CDSSNPS.txt
#166592

#Frequency calculations for each 
vcftools --vcf Ceratina.vcf --keep samples.txt --positions Matches.txt --freq --out TargetFreq
vcftools --vcf Ceratina.vcf --remove samples.txt --positions Matches.txt --freq --out OutgroupFreq

paste -d '\t' TargetFreq.frq OutgroupFreq.frq | sed 's/:/\t/g' |cut -f 1,2,6,8,14,16 | sed '1d' | sed '1iChrom\tPos\tTarget_Ref\tTarget_Alt\tOut_Ref\tOut_Alt' > MergedFreq.txt

awk 'BEGIN {IFS = OFS = "\t"} {if ($3 == 1 && $6 == 1)print $0, "FIXED"; 
else if ($4 == 1 && $5 == 1) print $0, "FIXED";
else if ($3 == 1 && $5 == 1) print $0, "SAME";
else if ($3 > 0 && $3 < 1 || $5 > 0 && $5 < 1) print $0, "POLY"; 
else print $0, "SAME"}'  MergedFreq.txt > LabeledVariants.txt

#Only interested in the ones in genes, and those that are fixed or polymorphic 
sed 's/\t/_/' LabeledVariants.txt | sed '/SAME/ d' > FixedPoly.txt
cut -f 2,3 CDSSNPS.txt | sed 's/\t/_/' > GenicLoci.txt

grep -wFf GenicLoci.txt FixedPoly.txt > GenicPolyFixed.txt

#SNP ids
awk 'BEGIN {IFS = OFS = "\t"} {print $2,$3,$1}' CDSSNPS.txt | sed 's/\t/_/' > SNPIds.txt

#Merge and aggregate in R
R
setwd('/data3/african/sideproject/dova')
id = read.table(file='SNPIds.txt', header=F)
names(id) = c("CHROM_POS", "GENE")
label = read.table(file="GenicPolyFixed.txt", header=F)
names(label) = c("CHROM_POS", "Target_Ref", "Target_Alt", "Out_Ref", "Out_Alt", "Label")

m = merge(label, id, by="CHROM_POS")
agg = aggregate(Label ~ GENE, data=m, summary)

write.table(m, "Categories_Gene.txt", sep='\t', quote=F, row.names=F)
write.table(agg, "SummaryTable.txt", sep='\t', quote=F, row.names=F)

#categories genes -- get LOC names and aggregate any overlapping genes - change any fixed/poly combinations to poly 

#in python -- determine synonymous and nonsynonymous 