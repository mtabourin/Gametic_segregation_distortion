# Search for male gametic segregation in F1 hybrids

This pipeline allows to search for male segregation distortions in F1 hybrids. It works in 4 scripts.

----------------

## Data

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



### Script_1.sh


### Script_2.sh


### Script_3.sh


### Script_4.sh

----------------

## Test
