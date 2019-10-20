**GALAXY FETCH SEQUENCES**

*gfetchseq [interval_file] [genome] [api_key] [upper? yes/no] [purge Galaxy? yes/no]*

input a gatk_interval file (.interval), a valid genome code ("hg19", "hg38", ...), your Galaxy API key, 
a flag to instruct the program to convert DNA alphabet to uppercase (yes/no),
and another flag to remove files from Galaxy (yes/no).

The output will include a new bed file converted from the interval source file, and a formatted fasta file containing all sequences, downloaded to current folder, with record IDs in the following format:

*genName_genStart-genEnd_repNameStrand*

DNA alphabet will be converted to upper case if the flag was "yes".

Condition: the initial .interval must be in one of these formets:

1) genName genStart genEnd
2) genName genStart genEnd Strand
3) genName genStart genEnd repName
4) genName genStart genEnd Strand repName
5) genName genStart genEnd repName Strand

|

**INSTALL**:

$ git clone http://github.com/lexxxxxxa/gfetchseq

**RUN**:

A) As command:

A-1) from source directory:

$ *./gfetchseq [interval_file] [genome] [api_key] [upper yes/no] [purge yes/no]*

A-2) as command from anywhere:

first copy executable from source directory to /usr/bin:

$ cp gfetchseq /usr/bin

then command can be used from anywhere:

$ *gfetchseq [interval_file] [genome] [api_key] [upper? yes/no] [purge? yes/no]*

|

B) As function:

place the *gfetchseq.py* in a folder along with an empty *__init__.py* file
import the module as *from folder.gfetchseq import gfetchseq*
then use the function in your script with *gfetchseq(interval_file, genome, api_key, yes/no, yes/no)*

|

**REQUIREMENTS**:

python >= 3.6.7

bioblend >= 0.13.0

biopython >= 1.74

pandas >= 0.25.1

termcolor >= 1.1.0
