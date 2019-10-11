**GALAXY FETCH SEQUENCES**

*gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no]*

input a gatk_interval file (.interval), a valid genome code ("hg19", "hg38", ...), your Galaxy API key, 
a flag to instruct the program to convert DNA alphabet to lowercase (yes/no),
and another flag to remove files from Galaxy (yes/no).

The output will be a friendly formatted fasta, with record IDs in the following form:

*genName_genStart-genEndStrand*

and adding in each record description the name (if applicable) and/or index of the sequence.

DNA alphabet converted to lower case if the flag was "yes".

|

**INSTALL**:

git clone http://github.com/lexxxxxxa/gfetchseq

**RUN**:

1) from source directory:

$ *./gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no]*

2) as command line:

add source dir to PATH and run from anywhere as:

$ *gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no]*

|

**Requirements**:

python >= 3.6.7

bioblend >= 0.13.0

biopython >= 1.74

pandas >= 0.25.1

termcolor >= 1.1.0
