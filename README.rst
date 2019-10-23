**GALAXY FETCH SEQUENCES**

*gfetchseq [interval/bed file] [genome] [api_key] [upper? yes/no] [purge Galaxy? yes/no]*

|

**INPUT**

1) Source file
 
1a) INTERVAL file (.interval) in one of the following formats:

genName genStart genEnd

genName genStart genEnd Strand

genName genStart genEnd repName

genName genStart genEnd Strand repName

genName genStart genEnd repName Strand

|

1b) BED file (.bed), in this format only: 

genName genStart genEnd repName Score Strand

|

2) a valid genome code ("hg19", "hg38", ...)

3) your Galaxy API key

4) a flag to instruct the program to convert DNA alphabet to uppercase (yes/no)

5) another flag to remove files from Galaxy (yes/no).

|

**OUTPUT**

The output will include a new bed file converted if the input file was interval, and a formatted fasta file containing all sequences, downloaded to current folder, with record IDs in the following format:

*genName_genStart-genEnd_repNameStrand*

DNA alphabet will be converted to upper case if the flag was "yes".


|

**INSTALL**:

$ git clone http://github.com/lexxxxxxa/gfetchseq

**RUN**:

A) As command:

copy executable from command directory to /usr/bin:

$ cp gfetchseq /usr/bin

use from anywhere:

$ *gfetchseq [interval/bed] [genome] [api_key] [upper? yes/no] [purge? yes/no]*

|

B) As function:

- create a folder at the root of your project, move the *gfetchseq.py* from the function directory, along with an empty *__init__.py* file
- import the module as *from folder.gfetchseq import gfetchseq*
- then use in script with *gfetchseq(interval/bed, genome, api_key, yes/no, yes/no)*

|

**REQUIREMENTS**:

python >= 3.6.7

bioblend >= 0.13.0

biopython >= 1.74

pandas >= 0.25.1

termcolor >= 1.1.0
