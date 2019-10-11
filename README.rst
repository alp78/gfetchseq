Galaxy Fetch Sequence

gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no] 

input a gatk_interval file (.interval), a valid genome code, your Galaxy API key, 
a flag to instruct the program to convert DNA alphabet to lowercase or not (yes/no),
and another flag to clean your Galaxy history after completion (yes/no).

The output will be a friendly formatted fasta, with record id [genName_genStart-genEndStrand],
and DNA alphabet converted to lower case if the flag was "yes".

install:

git clone http://github.com/lexxxxxxa/gfetchseq

run:

1) from source directory:
$ ./gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no]

2) as command line:
add source dir to PATH and run from anywhere as:

gfetchseq [interval_file] [genome] [api_key] [lower yes/no] [purge yes/no]



Requirements:

python >= 3.6.7

bioblend >= 0.13.0

biopython >= 1.74

pandas >= 0.25.1

termcolor >= 1.1.0
