#!/bin/bash
#
#SBATCH --ntasks-per-node=8
#SBATCH --partition=long
#SBATCH -n 16 #threads
#SBATCH --mem 80GB

module load samtools/1.9
module load parallel/20190322
module load gatk4/4.2.3.0
module load picard/2.23.5

#BASH SCRIPT FOR VARIANT CALLING


###Define variables
DIR=/shared/projects/gametic_segregation_distortion
DIRDATA="$DIR"/data
DIRINFO="$DIRDATA"/genealogies_table_TRUE.txt
DIRREF="$DIRDATA"/reference
REF=Alyrata_107

DIRRESULT="$DIR"/results
DIRMAP="$DIRRESULT"/mapping
DIRVCF="$DIRRESULT"/VCF
DIRVCFCROSS="$DIRRESULT"/VCF_cross
DIRTABLE="$DIRRESULT"/table

declare -a indivs
indivs=(
	'P01'
	'P02'
	'L'
	'P'
)

###Define the name of individus

##P01
name_P01=$(grep $1 "$DIRINFO" | grep -w P01)
name_P01=(${name_P01// / })

##P02
name_P02=$(grep $1 "$DIRINFO" | grep -w P02)
name_P02=(${name_P02// / })

##L
name_L=$(grep $1 "$DIRINFO" | grep -w L)
name_L=(${name_L// / })

##P
name_P=$(grep $1 "$DIRINFO" | grep -w p)
name_P=(${name_P// / })

declare -A names
names=( ['P01']="${name_P01[2]}" ['P02']="${name_P02[2]}" ['L']="${name_L[2]}" ['P']="${name_P[2]}" )


### Retrieval of chromosomes names
fasta_seq_ID=$(sed -n '/^>/p' "$DIRREF"/"$REF".fa | sed -e "s/^>//")
fasta_seq_ID=(${fasta_seq_ID//"\n"/ })

for indiv in ${indivs[@]}
do
	all_gz_picard=()
	all_gz=()
	all_tbi=()
	
	echo Variant Calling
	
	for chrom in "${fasta_seq_ID[@]}"
	do
	
		if [ ! -d "$DIRVCF"/$1 ]
		then 
			mkdir "$DIRVCF"/$1
		fi
	
		if [ ! -d "$DIRVCF"/$1/"$indiv" ]
		then 
			mkdir "$DIRVCF"/$1/"$indiv"
		fi
	
		all_gz_picard[${#all_gz_picard[@]}]="-I "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}."$chrom".sorted.MD.RG.g.vcf.gz"
		all_gz[${#all_gz[@]}]="$DIRVCF"/$1/"$indiv"/${names["$indiv"]}."$chrom".sorted.MD.RG.g.vcf.gz
		all_tbi[${#all_tbi[@]}]="$DIRVCF"/$1/"$indiv"/${names["$indiv"]}."$chrom".sorted.MD.RG.g.vcf.gz.tbi
	
	done
	
	(for i in "${fasta_seq_ID[@]}";do echo "$i";done) | parallel -j 8 "gatk HaplotypeCaller -R "$DIRREF"/"$REF".fa -I "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.RG.bam -O "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.{}.sorted.MD.RG.g.vcf.gz -ERC GVCF --native-pair-hmm-threads 2 -L {}"
	
	wait;

	echo Merge gVCF "$indiv"
	
	picard GatherVcfs ${all_gz_picard[@]} -O "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf
	
	bgzip -c "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf > "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf.gz
	tabix -p vcf "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf.gz
	
	if [ -s "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf.gz ]
	then 
		rm ${all_gz[@]} ${all_tbi[@]}
		rm "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf
	fi
done

### Create a database for the cross, with the gVCF, if it does not exist and add the gVCF if the database for the cross already exists

## Index the reference file
if [ ! -f "$DIRREF"/"$REF".bed ]
then 
	echo Index the reference file
	
	samtools faidx "$DIRREF"/"$REF".fa
	awk 'BEGIN {FS="\t"}; {print $1 FS "0" FS $2}' "$DIRREF"/"$REF".fa.fai > "$DIRREF"/"$REF".bed
		
fi

echo Create a database dor the cross "$1"

declare -a indivs
indivs=(
        'P01'
	'P02'
	'L'
	'P'
)

format4variant=()
for indiv in ${indivs[@]}
do
	format4variant[${#format4variant[@]}]="--variant "$DIRVCF"/$1/"$indiv"/${names["$indiv"]}.g.vcf.gz"
done

gatk GenomicsDBImport -R "$DIRREF"/"$REF".fa ${format4variant[@]} --genomicsdb-workspace-path "$DIRVCFCROSS"/$1.GATK.database --intervals "$DIRREF"/"$REF".bed

cd "$DIRVCFCROSS"

gatk GenotypeGVCFs -R "$DIRREF"/"$REF".fa -V gendb://$1.GATK.database -O "$DIRVCFCROSS"/$1.vcf.gz

echo VCF to Table

gatk VariantsToTable -R "$DIRREF"/"$REF".fa -V "$DIRVCFCROSS"/$1.vcf.gz -F CHROM -F POS -F TYPE -GF GT -GF GQ -GF DP -GF AD -O "$DIRTABLE"/$1.table










