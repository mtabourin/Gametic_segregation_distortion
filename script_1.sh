#!/bin/bash
#
#SBATCH --cpus-per-task=45
#SBATCH --partition=long

module load trimmomatic/0.39
module load bwa/ #Version 0.7.17-r1188
module load samtools/1.9 
module load picard/2.23.5

#BASH SCRIPT FOR PREPARING DATA FOR VARIANT CALLING


###Define variables
DIR=/shared/projects/gametic_segregation_distortion
DIRDATA="$DIR"/data
DIRINFO="$DIRDATA"/genealogies_table_TRUE.txt
DIRREF="$DIRDATA"/reference
REF=Alyrata_107

DIRRESULT="$DIR"/results
DIRTMP="$DIRRESULT"/tmp
DIRTRIM="$DIRRESULT"/trimmed_reads
DIRMAP="$DIRRESULT"/mapping

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


for indiv in ${indivs[@]}
do
	##Pooling of fastq files if necessary
	
	if (( $(ls "$DIRDATA"/$1/"$indiv" | wc -l) > 2 ))
	then
		echo Pooling fastq files for "$indiv"
		
		##File decompression
		
		for file in $(ls "$DIRDATA"/$1/"$indiv")
		do
			name=$(basename "$file" .gz)
			pigz -d --force -c -p $SLURM_CPUS_PER_TASK "$DIRDATA"/$1/"$indiv"/"$file" > "$DIRTMP"/"$name"
		done
		
		##Pooling of R1 fastqc file
		
		file_R1=$(ls "$DIRDATA"/$1/"$indiv" | grep -i "R1")
		file_R1=(${file_R1//"\n"/ })
		p1=""
		
		for ((i=0;i<=${#file_R1[@]}-1;i++))
		do
			p1+="$DIRTMP"/$(basename "${file_R1[$i]}" .gz)
			p1+=" "
		done
		name=${file_R1[0]%.*}
		cat $p1 > $DIRDATA/$1/$indiv/${file_R1[0]%.*}_pooled.fastq
		pigz -p $SLURM_CPUS_PER_TASK $DIRDATA/$1/$indiv/${file_R1[0]%.*}_pooled.fastq
		
		##Pooling of R2 fastqc file
		
		file_R2=$(ls "$DIRDATA"/$1/"$indiv" | grep -i "R2")
		file_R2=(${file_R2//"\n"/ })
		p2=""
		
		for ((i=0;i<=${#file_R2[@]}-1;i++))
		do
			p2+="$DIRTMP"/$(basename "${file_R2[$i]}" .gz)
			p2+=" "
		done
		name=${file_R2[0]%.*}
		cat $p2 > $DIRDATA/$1/$indiv/${file_R2[0]%.*}_pooled.fastq
		pigz -p $SLURM_CPUS_PER_TASK $DIRDATA/$1/$indiv/${file_R2[0]%.*}_pooled.fastq
		
		##Delete temporary files
		
		rm "$DIRTMP"/*	
	fi
	
	##Trimming of fastq files
	
	echo Trimming fastq files for "$indiv"
	
	if [ ! -d "$DIRTRIM"/$1 ]
	then 
		mkdir "$DIRTRIM"/$1
	fi
	
	if [ ! -d "$DIRTRIM"/$1/"$indiv" ]
	then 
		mkdir "$DIRTRIM"/$1/"$indiv"
	fi
	
	if (( $(ls "$DIRDATA"/$1/"$indiv" | wc -l) > 2 ))
	then
		file_pooled=$(ls "$DIRDATA"/$1/"$indiv" | grep -i "pooled")
		file_pooled=(${file_pooled//"\n"/ })
		
		trimmomatic PE -threads $SLURM_CPUS_PER_TASK "$DIRDATA"/$1/"$indiv"/${file_pooled[0]} "$DIRDATA"/$1/"$indiv"/${file_pooled[1]} "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1_unpaired.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2_unpaired.fastq.gz ILLUMINACLIP:"$DIRDATA"/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75
		
		rm -f "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1_unpaired.fastq.gz
		rm -f "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2_unpaired.fastq.gz
	
	else
		file=$(ls "$DIRDATA"/$1/"$indiv")
		file=(${file//"\n"/ })
		
		trimmomatic PE -threads $SLURM_CPUS_PER_TASK "$DIRDATA"/$1/"$indiv"/${file[0]} "$DIRDATA"/$1/"$indiv"/${file[1]} "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1_unpaired.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2_unpaired.fastq.gz ILLUMINACLIP:"$DIRDATA"/adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75 
		
		rm -f "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1_unpaired.fastq.gz
		rm -f "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2_unpaired.fastq.gz
		
	fi
	
	##Mapping
		
	##Index the reference file
	if [ ! -f "$DIRREF"/"$REF".fa.bwt ]
	then 
		echo Index the reference file
		
		bwa index "$DIRREF"/"$REF".fa
		samtools faidx "$DIRREF"/"$REF".fa
		picard CreateSequenceDictionary -R "$DIRREF"/"$REF".fa -O "$DIRREF"/"$REF".dict
	fi
			
	##Do mapping
		
	echo Mapping for "$indiv"
	
	if [ ! -d "$DIRMAP"/$1 ]
	then 
		mkdir "$DIRMAP"/$1
	fi
	
	if [ ! -d "$DIRMAP"/$1/"$indiv" ]
	then 
		mkdir "$DIRMAP"/$1/"$indiv"
	fi	
	
	bwa mem -t $SLURM_CPUS_PER_TASK "$DIRREF"/"$REF".fa "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R1.fastq.gz "$DIRTRIM"/$1/"$indiv"/${names["$indiv"]}_R2.fastq.gz > "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sam
				
	##Change file format
		
	echo Change file format for "$indiv"
	
	file=$(ls "$DIRDATA"/$1/"$indiv")
	file=(${file//"\n"/ })
		
	name_PU=$(less "$DIRDATA"/$1/"$indiv"/${file[0]} | head -n 1)
	name_PU=(${name_PU//:/ })
	name_PU=${name_PU[2]}${name_PU[3]}

	
	samtools view -@ $SLURM_CPUS_PER_TASK -Sb -h "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sam > "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.bam
	
	samtools sort -@ $SLURM_CPUS_PER_TASK "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.bam > "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.bam
	
	java -XX:ParallelGCThreads=$SLURM_CPUS_PER_TASK -jar /shared/ifbstor1/software/miniconda/envs/picard-2.23.5/share/picard-2.23.5-0/picard.jar MarkDuplicates -I "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.bam -O "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.bam -M trash
	
	java -XX:ParallelGCThreads=$SLURM_CPUS_PER_TASK -jar /shared/ifbstor1/software/miniconda/envs/picard-2.23.5/share/picard-2.23.5-0/picard.jar AddOrReplaceReadGroups -I "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.bam -O "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.RG.bam -RGID ${names["$indiv"]} -RGLB "$indiv" -RGPL ILLUMINA -RGPU "$name_PU" -RGSM ${names["$indiv"]} 
	
	rm -f "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sam
	rm -f "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.bam
	rm -f "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.bam
	rm -f "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.bam
	
	samtools index "$DIRMAP"/$1/"$indiv"/${names["$indiv"]}.sorted.MD.RG.bam
done	
	

