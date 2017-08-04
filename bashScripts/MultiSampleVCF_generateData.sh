#!/bin/bash
#
#PBS -r n

##Job settings
#PBS -N 20_variantList
#PBS -m be
#PBS -q default
#PBS -o Logs/20_variantList_stdout.$PBS_JOBID
#PBS -e Logs/20_variantList_stderr.$PBS_JOBID
###PBS -M r.govind@rbht.nhs.uk

##Job configuration
#PBS -V

##Job Resources
#PBS -l nodes=1:ppn=56
#PBS -l mem=20gb
#PBS -l walltime=120:00:00

#select gene_symbol,gene_ensg from ICC_MUTATIONS.gene where gene_symbol in ('BAG3','LMNA','TCAP','TNNC1','TNNT2','TPM1','DSP','SCN5A','TTN','VCL','MYH7');

source ./CONFIG

group=$TopDir/ALLSamples

$Java7 $memory -jar $GATKs/GenomeAnalysisTK.jar \
		-R $Ref \
		-T VariantFiltration \
		-o $group.QDFilter.vcf \
		--variant $inFile \
		--filterExpression "QD < 4.0" \
		--filterName QDFilter \

	
#Remove filtered variants - grep PASS 

(grep ^"#" $group.QDFilter.vcf; grep PASS $group.QDFilter.vcf) > $group.PassOnly.vcf 
	
# TITIN PSI - constititive exons 
TTN_constitutive_exons=$TopDir/TTN_constitutive_exons.bed

#sort VCF Files
perl $vcftools/vcf-sort -c $group.PassOnly.vcf  > $group.PassOnly.sorted.vcf

#zip VCF File
$samtools/bgzip -c $group.PassOnly.sorted.vcf > $group.PassOnly.sorted.vcf.gz

#Create tbix index
$samtools/tabix -p vcf $group.PassOnly.sorted.vcf.gz

#tabix
$samtools/tabix $group.PassOnly.sorted.vcf.gz -R $TTN_constitutive_exons > $group.PassOnly.TTN_constitutive_exons.vcf.gz


# Unzip TTN_constitutive_variants, get non-chr2 variants, combine two
$samtools/bgzip -d $group.PassOnly.TTN_constitutive_exons.vcf.gz
(grep ^"#" $group.PassOnly.vcf; grep -v chr2 $group.PassOnly.vcf) > $group.PassOnly.Chr2_Removed.vcf 
#grep -wv ^"chr2" $group.PassOnly.vcf > $group.PassOnly.Chr2_Removed.vcf 


cat $group.PassOnly.Chr2_Removed.vcf $group.PassOnly.TTN_constitutive_exons.vcf > $group.PassOnly.TTN_constitutive_exons.AllChr.vcf 
perl $vcftools/vcf-sort -c $group.PassOnly.TTN_constitutive_exons.AllChr.vcf   > $group.PassOnly.TTN_constitutive_exons.AllChr.sorted.vcf 

Filtered=$group.PassOnly.TTN_constitutive_exons.AllChr.sorted.vcf

#Tablize
		
python $LOFTEE/tableize_vcf.py \
--vcf $Filtered --out $out_tableize --split_by_transcript --all_csqs --do_not_minrep --include_id --info ABHet,ABHom,AC,AF,AN,DP,FS,HaplotypeScore,MLEAC,MLEAF,MQ,MQ0,QD \
--vep_info Consequence,SYMBOL,SYMBOL_SOURCE,Gene,Feature,BIOTYPE,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variation,STRAND,CANONICAL,CCDS,ENSP,SIFT,PolyPhen,ExAC_MAF,PUBMED \
--samples --mysql


# Filter by genes on Interest 

#SELECT gene_symbol,gene_canon_transcript FROM ICC_MUTATIONS.gene WHERE gene_symbol IN ("BAG3","LMNA","TCAP","TNNC1","TNNT2","TPM1","DSP","SCN5A","TTN","VCL","MYH7");
TranscriptList="CHROM|ENST00000333535|ENST00000232975|ENST00000368300|ENST00000211998|ENST00000355349|ENST00000367318|ENST00000403994|ENST00000309889|ENST00000589042|ENST00000379802|ENST00000369085"

grep -E $TranscriptList $out_tableize > $TranscriptsOfInterest

#Filter by VarType : Only look at Truncating / Non Truncating Variants
# awk '{print $21}' $TopDir/ALLSamples.TranscriptsOfInterest.txt | sed 's/&/\n/g' | sort | uniq -c | sort -nr
   # 1201 missense_variant            : keep : non-truncating 
    # 596 synonymous_variant          : discard
    # 302 intron_variant              : discard 
     # 70 splice_region_variant       : [this annotation comes with oother annotations like missense_variant,intron_variant,etc ]
     # 56 frameshift_variant          : keep : truncating 
     # 49 stop_gained                 : keep : truncating 
     # 34 3_prime_UTR_variant         : discard
     # 11 splice_donor_variant        : keep : truncating
      # 9 inframe_deletion            : keep : non-truncating
      # 8 5_prime_UTR_variant         : discard
      # 7 downstream_gene_variant     : discard
      # 4 splice_acceptor_variant     : keep : truncating
      # 3 inframe_insertion           : keep : non-truncating
      # 1 protein_altering_variant    : keep : non-truncating
      # 1 coding_sequence_variant     : [this annotation comes with oother annotations like splice_donor_variant,splice_acceptor_variant, etc ]

# GENE	Non-Truncating	Truncating	Transcript ID
# BAG3	YES	            YES	         ENST00000369085
# LMNA	YES	            Yes	         ENST00000368300
# TCAP	YES	            Yes	         ENST00000309889
# TNNC1	YES	            Yes	         ENST00000232975
# TNNT2	YES	            Yes	         ENST00000367318
# TPM1	YES	            Yes	         ENST00000403994
# DSP	NO	            Yes	         ENST00000379802
# SCN5A	NO	            Yes	         ENST00000333535
# TTN	NO	            Yes	         ENST00000589042
# VCL	NO	            Yes	         ENST00000211998
# MYH7	Yes	            NO	         ENST00000355349

	  
PAV_List="CHROM|missense_variant|frameshift_variant|stop_gained|inframe_deletion|splice_donor_variant|splice_acceptor_variant|inframe_insertion|protein_altering_variant"

TV_List="CHROM|frameshift_variant|stop_gained|splice_donor_variant|splice_acceptor_variant"

nTV_List="CHROM|missense_variant|inframe_deletion|inframe_insertion|protein_altering_variant"
	  
(grep BAG3 $TranscriptsOfInterest | grep -E $PAV_List; grep LMNA $TranscriptsOfInterest | grep -E $PAV_List; grep TCAP $TranscriptsOfInterest | grep -E $PAV_List; grep TNNC1 $TranscriptsOfInterest | grep -E $PAV_List; grep TNNT2 $TranscriptsOfInterest | grep -E $PAV_List; grep TPM1 $TranscriptsOfInterest | grep -E $PAV_List; grep DSP $TranscriptsOfInterest | grep -E $TV_List; grep SCN5A $TranscriptsOfInterest | grep -E $TV_List; grep TTN $TranscriptsOfInterest | grep -E $TV_List; grep VCL $TranscriptsOfInterest | grep -E $TV_List; grep MYH7 $TranscriptsOfInterest | grep -E $nTV_List; ) > $PAV
	  

# awk to form smaller table 

awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$21"\t"$22"\t"$25"\t"$27"\t"$28"\t"$20}' $PAV | sort -k6 > $TopDir/ALLSamples_Table.txt

# Add ExAC and FAF frequencies 
#get subset of ExAC [only genes on interest] #first get regions we need in bed File


grep -E "BAG3|LMNA|TCAP|TNNC1|TNNT2|TPM1|DSP|SCN5A|TTN|VCL|MYH7" $ICC_169Genes | sed 's/chr//g' > $TopDir/GenesOfInterest.bed

$samtools/tabix $ExAC -R $TopDir/GenesOfInterest.bed > $TopDir/ExAC.r0.3.1.GenesOfInterest.vcf.gz

ExAC=$TopDir/ExAC.r0.3.1.GenesOfInterest.vcf.gz

#zgrep -v '^#' $ExAC | awk '{ if ($1=="10" && $2=="121436286" && $4=="C" && $5=="T")  print $1";"$2";"$4";"$5";"$8 }' | awk -F";" '{for(i=1;i<=NF;i++){if ($i ~ /^AF=/){print $1"\t"$2"\t"$3"\t"$4"\t"$i}}}'

echo CHROM POS REF ALT CSQ GENE TRANSCRIPT HGVSc HGVSp SampleID ExAC_AF af_filter af_filter_pop DCM_max_allele_freq > $TopDir/ALLSamples_Table.freq.txt

#From disease table in ICC_MUTATIONS
DCM_max_allele_freq=0.00008400

while read line; do

	CHROM=`echo $line | awk '{print $1}'`
	POS=`echo $line | awk '{print $2}'`
	REF=`echo $line | awk '{print $3}'`
	ALT=`echo $line | awk '{print $4}'`
	
	CHROM=`echo $CHROM | sed 's/chr//' `

	ExAC_AF=`zgrep -v '^#' $ExAC | awk -v C="$CHROM" -v P="$POS" -v R="$REF" -v A="$ALT"  '{ if ($1==C && $2==P && $4==R && $5==A )  print $1";"$2";"$4";"$5";"$8 }' | awk -F";" '{for(i=1;i<=NF;i++){if ($i ~ /^AF=/){print $i}}}' | sed 's/AF=//'`

	af_filter_AF=`zcat $af_filter | awk -v C="$CHROM" -v P="$POS" -v R="$REF" -v A="$ALT"  '{ if ($1==C && $2==P && $3==R && $4==A )  print $5"\t"$6 }' `
		
	if [[ $ExAC_AF == "" ]]; then
		ExAC_AF_temp=`zgrep -v '^#' $ExAC | awk -v C="$CHROM" -v P="$POS" -v R="$REF" -v A="$ALT" '{ if ($1==C && $2==P && $4==R )  print $1";"$2";"$4";"$5";"$8 }' | awk -F";" '{for(i=1;i<=NF;i++){if ($i ~ /^AF=/){print $i}}}' | sed 's/AF=//'`
		ExAC_ALT=`zgrep -v '^#' $ExAC | awk -v C="$CHROM" -v P="$POS" -v R="$REF" -v A="$ALT" '{ if ($1==C && $2==P && $4==R )  print $5 }'`	
		
		if [[ $ExAC_AF_temp != "" ]]
		then
			IFS=', ' read -r -a alt_array <<< "$ExAC_ALT"
			IFS=', ' read -r -a af_array <<< "$ExAC_AF_temp"
			
			for index in "${!alt_array[@]}"
			do
				alt_temp=${alt_array[index]}
				af_temp=${af_array[index]}
				
				if [[ $alt_temp == $ALT ]]; then
					ExAC_AF=$af_temp
				fi
				
			done			
		fi
	fi
	
	if [[ $ExAC_AF == "" ]]; then
		ExAC_AF=0
	fi	
	
	if [[ $af_filter_AF == "" ]]; then
		af_filter_AF="0\t-"
	fi		
	
	echo -e "$line $ExAC_AF $af_filter_AF $DCM_max_allele_freq" >> $TopDir/ALLSamples_Table.freq.txt

done < $TopDir/ALLSamples_Table.txt  

perl -p -i -e 's/ /\t/g' $TopDir/ALLSamples_Table.freq.txt


# make file per Sample 

while read line; do
	
	variant=`echo $line | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$11"\t"$12"\t"$13"\t"$14}'`
	samples=`echo $line | awk '{print $10}' | sed 's/,/ /g'`
	
	for sample in $samples
	do
		if [[ $samples = "SampleID" ]]; then
			echo CHROM POS REF ALT CSQ GENE TRANSCRIPT HGVSc HGVSp ExAC_AF af_filter af_filter_pop DCM_max_allele_freq SampleID > $TopDir/ALLSamples_Table.perSample.txt
		else		
			echo -e "$variant\t$sample" >> $TopDir/ALLSamples_Table.perSample.txt
		fi
	done

done < $TopDir/ALLSamples_Table.freq.txt

perl -p -i -e 's/ /\t/g' $TopDir/ALLSamples_Table.perSample.txt


# add VarType Column - ie TV or nTV

while read line; do
	
	if echo $line | grep -qE CHROM ; then
		echo $line Var_Group > $TopDir/ALLSamples_Table.VarGroup.txt
		
	elif echo $line | grep -qE $TV_List; then
		echo $line TV >> $TopDir/ALLSamples_Table.VarGroup.txt
		
	elif echo $line | grep -qE $nTV_List; then
		echo $line non_TV >> $TopDir/ALLSamples_Table.VarGroup.txt
		
	else
		echo $line unknown >> $TopDir/ALLSamples_Table.VarGroup.txt
	fi

done < $TopDir/ALLSamples_Table.perSample.txt

perl -p -i -e 's/ /\t/g' $TopDir/ALLSamples_Table.VarGroup.txt

# add Cohort Group

while read line; do
	sample=`echo $line | awk '{print $14}'`
	
	if [[ $sample == SampleID ]]; then
		echo $line group > $TopDir/ALLSamples_Table.group.txt
	elif [[ $sample == 10* ]]; then
		echo $line DCM >> $TopDir/ALLSamples_Table.group.txt
	elif [[ $sample == 11* ]]; then
		echo $line DCM >> $TopDir/ALLSamples_Table.group.txt
	elif [[ $sample == 12* ]]; then
		echo $line DCM >> $TopDir/ALLSamples_Table.group.txt
	elif [[ $sample == 14* ]]; then
		echo $line HVOL >> $TopDir/ALLSamples_Table.group.txt
	elif [[ $sample == 20* ]]; then
		echo $line ACM >> $TopDir/ALLSamples_Table.group.txt
	else
		echo $line unknown >> $TopDir/ALLSamples_Table.group.txt		
	fi
	
done < $TopDir/ALLSamples_Table.VarGroup.txt

perl -p -i -e 's/ /\t/g' $TopDir/ALLSamples_Table.group.txt

# Remove Filtered Out Samples  

Samples2Keep=`awk '{print $1}' $Filtered_SampleID | tr '\n' '|'`
Samples2Keep=$Samples2Keep"CHROM"

grep -E $Samples2Keep $TopDir/ALLSamples_Table.group.txt > $TopDir/FilteredSamples.variantList.txt

# Remove Common Variants 

		echo "Variants <- read.table(\"$TopDir/FilteredSamples.variantList.txt\", header=T)
			RareVar=Variants[which(Variants\$af_filter<Variants\$DCM_max_allele_freq),]
			write.table(RareVar, quote=FALSE, sep =\"\t\", row.names = F, file=\"$TopDir/Final_VariantList_Rare.txt\") 
			" > $TopDir/RemoveCommonVariants.r

		$RSoftware CMD BATCH $TopDir/RemoveCommonVariants.r
		
cp $TopDir/Final_VariantList_Rare.txt $Rproject_input_files_DIR/Variants_filteredbyGeneVarType_allGroups_allfreq.txt

		
# Seperate into groups 

#awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$21"\t"$22"\t"$25"\t"$27"\t"$28"\t"$20}' $PAV | grep -E $geneList | sort -k9 | awk '$9 ~ /^20/' | sort -k6 > $TopDir/ACM_variants.txt


