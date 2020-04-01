# fastv
fastv: an ultra-fast sequencing data analysis tool for screening of viral infectious diseases like COVID-19

# what's fastv?
`fastv` is designed to detect SARS-CoV-2 from sequencing data in FASTQ format, but will be also able to detect other viral sequences. `fastv` accepts of RAW FASTQ files (both SE/PE), performs quality filtering (cut adapters, remove low quality reads), and then scan for the viral sequences.

# take a quick glance of the informative report
* Sample HTML report: http://opengene.org/fastv/fastv.html
* Sample JSON report: http://opengene.org/fastv/fastv.json

# try fastv to generate above reports
* FASTQ file for testing: http://opengene.org/fastv/testdata.fq.gz
* Command for testing: 
```shell
./fastv -i testdata.fq.gz
``

# quick examples
Single-end data
```shell
fastv -i testdata.fq.gz
```
paired-end data
```shell
fastv -i R1.fq.gz -I R2.fq.gz
```
output the reads containing target viral sequences
```shell
fastv -i R1.fq.gz -I R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz
```
specify the KMER file and genome file
```shell
fastv -i R1.fq.gz -I R2.fq.gz -k data/SARS-CoV-2.kmer.fa -g data/SARS-CoV-2.genomes.fa
```

# get fastv
## download binary 
This binary is only for Linux systems: http://opengene.org/fastv/fastv
```shell
# this binary was compiled on CentOS, and tested on CentOS/Ubuntu
wget http://opengene.org/fastv/fastv
chmod a+x ./fastv
```
## or compile from source
```shell
git clone https://github.com/OpenGene/fastv.git

# step 3: build
cd fastv
make
```
