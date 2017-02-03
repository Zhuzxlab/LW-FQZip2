# LW-FQZip2

LW-FQZip 2 is a lossless reference-based compression method targeting FASTQ files. It improved from the light-weight reference based compression tool LW-FQZip 1 (Y. Zhang et al., BMC Bioinformatics, 16:188, 2015) by introducing more efficient coding schemes and parallelism. LW-FQZip 2 is capable of obtaining superior compression ratios at reasonable time cost and memory consumption. The competence enables LW-FQZip 2 to serve as a candidate tool for archival or storage space-sensitive applications of real-world sequencing data. 

### Download & Installation

LW-FQZip 2 is developed in C/C++ and runs on 32/64 bit Linux or Mac OS. A minimum memory of 1GB is suggested for good user experience. To install the program, please download and decompress the source code and then compile it in a GNU environment equipped with GCC compiler using the following commands:

	git clone https://github.com/Zhuzxlab/LW-FQZip2.git
	cd LW-FQZip2
	make clean
	make

The executable files ‘LWFQZip2’, ‘LWMapping’ and 'FQZip' will be generated in the same directory of the source code. The light-weight mapping model can be used independently with LW-FQZip 2. Make sure the archiving tool tar(http://www.gnu.org/software/tar/tar.html) have been installed and configured correctly before running LW-FQZip 2. High compression mode depends on compression algorithms: zpaq702(http://mattmahoney.net/dc/zpaq702.zip). To ensure that executable file zpaq has permission to perform.

### Commands

A simple example is provided as follows to illustrate the use of "LWFQZip2". To compress the sample FASTQ file SRR1063349.fastq, the command "LWFQZip2 -c" is executed with a reference NC_017634.1.fasta.

	LWFQZip2 -c -i SRR1063349.fastq -r NC_017634.1.fasta

where the target FASTQ file is first mapped to the reference obtaining an intermediate output file "SRR1063349.fastq.map.txt"(in SAM format) in the same directory, and then the original data is compressed based on this mapping results. A compressed file "SRR1063349.fastq.lz" is obtained.

To decompress the file, the command "LWFQZip2 -d" should be called.

	LWFQZip2 -d -i SRR1063349.fastq.lz -r NC_017634.1.fasta

More parameters can be specified for the mapping and compressiong parts as follows: 

COMMANDS AND OPTIONS
LWFQZip2 	<mode>...[options]
  	
    Mode:
  	
    -c 	compression.
    -d 	decompression.
    
    Compression/Decompression Options:
  	
    -i 	input FASTQ file or compressed file.
    -r 	input Reference file.
    -m 	maximal read length,ranging from 30000 to 300000 (Default: '-m 300000').
    -h 	help.
    -g 	high compression mode (slow), usage: LWFQZip2 -c -i input -r reference -g.
    -a 	assemble-based mode. An optional amount (Default: 0.3% of the original file size) of reads, which contain the predefined          
        prefix (Default: 'CG'), are assembled to form an artificial reference. At the end of the package, this artificial 
        reference is included.
        usage:
        LWFQZip2 –c –i input.fastq –a 0.003(default)
        LWFQZip2 –d –i input.fastq.lz -a.
    -s 	display the counts of the common prefixes in the reference, e.g. CG,AT,ATA...
    -v 	version.
  	
    Mapping Options:
  	
    -b 	the number of mapping thread(Default: 10, mininum: 1)
    -p 	specify the kmer prefixes, e.g.,'CG', 'AT', and 'TAG' (Default: '-p CG'). 'AA' is not recommended as a prefix.
    -k 	length of a kmer used in locate local alignment. (Default: '-k 8').
    -l 	the mini length of a legal alignment.(Default: '-l 12').
    -o 	the complementary palindrome mode.(Default: '-o 1' means start,otherwise'-o 0').
    -e 	the tolerance rate of mismatches.(Default: '-e 0.05').


More discussions are available in our paper and our page(http://csse.szu.edu.cn/staff/zhuzx/lwfqzip2/index.html). The details of implementation, data sets, and experimental studies are provided in the supplementary file. 
