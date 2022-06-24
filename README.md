# Search for male gametic segregation in F1 hybrids within the genus *Arabidopsis*

This pipeline allows to search for male segregation distortions in F1 hybrids. It works in 4 scripts.

----------------

## Data

The data used are available on the NCBI in the BioProject ID PRJNA750331.   
We used the reference genome of *Arabidopsis lyrata*[1].  

For each crossing, the data should be organized like this:

* data
  * name_cross  
    * P01  
    R1.fastq.gz   
    R2.fastq.gz   

    * P02  
    R1.fastq.gz   
    R2.fastq.gz   

    * L  
    R1.fastq.gz   
    R2.fastq.gz   

    * P   
    R1.fastq.gz  
    R2.fastq.gz   

With P01 for parent 1 data, P02 for parent 2 data, L for individual F1 somatic data, and P for individual F1 germline data.

A "genealogies_table_TRUE.txt" file is needed in the form:

| name_cross | sample | name_individual | name_data |
| ---------- | -------| --------------- | --------- |
| X01        | P01    | Pais09          | D111_EKDL200000635-1a-AK5849-AK5791_H2L72DSXY_L4 |
| X01        | P02    | Wall02          | 652_EKDL190132027-2a-AK6504-AK6658_HJFGTDSXX_L2  |
| X01        | L      | X01_L           | FPais09xWall024-2                                |
| X01        | P      | X01_p           | PPais09xWall024-2                                |
| X02        | P01    | Wall10          | D654_EKDL200000637-1a-AK6524-AK6699_H2KGYDSXY_L4 |
| ...        | ...    | ...             | ...                                              |

----------------

## Running

The scripts were coded for use with SLURM. The directory paths in each script must be modified, before using them, according to your folder hierarchy.  
The versions of the tools used are noted in the scripts.

### Script_1.sh

This first script cleans the raw data, aligns it to the reference genome and puts the alignment file in the right format to be able to perform variant analysis.  

To run:  
```
sbatch script_1.sh name_cross
```

### Script_2.sh

This script allows you to call variants to genotype single nucleotide polymorphisms and put the file in an easy-to-analyze format.

To run:    	
```
sbatch script_2.sh name_cross
```

### Script_3.sh

This script makes it possible to filter the table and to create files for the use of MAP_SD and the creation of graphs.

To run:     	
```
sbatch script_3.sh name_cross
```

### Script_4.sh

This script allows the creation of graphs before the use of MAP_SD and after the results of MAP_SD.

To run:	    
```
sbatch script_4.sh name_cross
```

----------------

## Results

In the table folder of the results folder, the filtered and unfiltered tables are available for our crosses.  

----------------

## References

[1] C. Roux, V. Castric, M. Pauwels, S. I. Wright, P. Saumitou-Laprade, et X. Vekemans, « Does Speciation
between Arabidopsis halleri and Arabidopsis lyrata Coincide with Major Changes in a Molecular Target of
Adaptation? », PLOS ONE, vol. 6, n o 11, p. e26872, nov. 2011, doi: 10.1371/journal.pone.0026872.
