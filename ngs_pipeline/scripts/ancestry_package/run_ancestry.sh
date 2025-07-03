#!/bin/bash

working_dir="/mnt/MoCha-NGS/active/changtn/PDX_ancestry"
data_dir="/mnt/projects/MoCha-hiseq-retro-scratch/BW_transfers/processedDATA"
out_dir="/mnt/MoCha-NGS/active/changtn/PDX_ancestry/output"
samplelist=$1


#sbatch -p norm -t 5:00:00 --nodes=2 --ntasks-per-node=16 --mem=10G run_ancestry.sh

#snpwt.bed => pdxnet_v2.snpwt_hg19.bed
#snpwt.NA => pdxnet_v2.snpwt_hg19_NA

module load samtools/1.8
module load bcftools/1.16

for samID in `cat $samplelist`;
do
	echo -e "Working on $samID"
	pat=`echo $samID|cut -f1 -d"~"`

	#bk
	samtools mpileup -q 30 -Q 20 -v -f $working_dir/ucsc.hg19.fasta -l $working_dir/pdxnet_v2.snpwt_hg19.bed $data_dir/$pat/20170910/$samID/${samID}.bwa.final.bam | bcftools call -c -Ov | bcftools filter -e 'ALT="."' > $out_dir/${samID}.ancestry_temp.vcf
       	bcftools annotate -c CHROM,FROM,TO,ID -a $working_dir/pdxnet_v2.snpwt_hg19.bed.gz $out_dir/${samID}.ancestry_temp.vcf > $out_dir/${samID}.ancestry.vcf

	#samtools mpileup -q 30 -Q 20 -f $working_dir/ucsc.hg19.fasta -l $working_dir/pdxnet_v2.snpwt_hg19.bed $data_dir/$pat/20170910/$samID/${samID}.bwa.final.bam | bcftools call -c -Ov | bcftools filter -e 'ALT="."'| bcftools annotate -c CHROM,FROM,TO,ID -a $working_dir/pdxnet_v2.snpwt_hg19.bed.gz > $out_dir/${samID}.ancestry.vcf

	python2.7 $working_dir/tools/gdc/vcf2eigenstrat.py -v $out_dir/${samID}.ancestry.vcf -o $out_dir/${samID}.ancestry

	echo "geno:   $out_dir/${samID}.ancestry.geno"     >$out_dir/${samID}.ancestry.tmp
	echo "snp:    $out_dir/${samID}.ancestry.snp"     >>$out_dir/${samID}.ancestry.tmp
	echo "ind:    $out_dir/${samID}.ancestry.ind"     >>$out_dir/${samID}.ancestry.tmp
	echo "snpwt:  $working_dir/pdxnet_v2.snpwt_hg19_NA"             >>$out_dir/${samID}.ancestry.tmp
	echo "predpcoutput:   $out_dir/${samID}.ancestry" >>$out_dir/${samID}.ancestry.tmp


	python2.7 $working_dir/inferancestry.py   --par $out_dir/${samID}.ancestry.tmp
	echo -e "#Sample ID\\tPopulation label\\tNumber of SNPs\\tPredicted PC1\\tPredicted PC2\\tPredicted PC3\\tPredicted PC4\\t% EUR ancestry\\t% AMR ancestry\\t% AFR ancestry\\t% SAS Ancestry\\t% EAS Ancestry" >>$out_dir/${samID}.ancestry
	sort $out_dir/${samID}.ancestry -o $out_dir/${samID}.ancestry
	rm $out_dir/${sample_ID}.ancestry.tmp $out_dir/${sample_ID}.ancestry.geno $out_dir/${sample_ID}.ancestry.snp $out_dir/${sample_ID}.ancestry.ind
	echo -e "Done for $samID"
done
