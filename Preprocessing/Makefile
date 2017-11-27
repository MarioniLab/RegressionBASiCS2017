#
# Recipes for the result files of the scRNAseq analysis of CD4 cells of mice
#


#### Applications used
mapper = /nfs/research2/marioni/Nils/software/gsnap/bin/bin/
bsub = Tdist/Scripts/bsub -K
samtools = /nfs/research2/marioni/Nils/software/samtools-1.2
htseq_count = /nfs/research2/marioni/Nils/software/HTSeq-0.6.1p1/scripts

#### Other variables
memlimit = 32000


# Merged files
merge = /nfs/research2/marioni/Nils/scImmu/data/merge

# Input
cur_input = /nfs/research2/marioni/Nils/BASiCS/Tcell_data/fastq_files/PRJEB12793/

# Output
cur_results = /nfs/research2/marioni/Nils/BASiCS/Tcell_data/bam_files

# Results
results = /nfs/research2/marioni/Nils/BASiCS/Tcell_data/results

#### Index
gsnap_index = /nfs/research2/marioni/Nils/software/gsnap/database

#### Mapping
mapped_reads = $(addprefix ${cur_results}/,$(patsubst %_1.fastq.gz,%.bam,$(notdir ${fastq_files})))
fastq_files = $(shell cat all-files-fastq.txt)

#### Counting
read_count = $(addprefix ${results}/,$(patsubst %.bam,%.txt,$(notdir ${bam_files})))
bam_files = $(shell cat all-files-bam.txt)

#### Mapping
# Generate File for files to map
all-files-fastq.txt: 
	ls ${cur_input}*/*1.fastq.gz > $@

all-files-bam.txt:
	ls ${cur_results}/*.bam > $@

.PHONY: mapped-reads
mapped-reads: ${mapped_reads}

# For the mapping adjust the reference!
${cur_results}/%.bam: ${cur_input}/*/%_1.fastq.gz ${cur_input}/*/%_2.fastq.gz  
	${bsub} -M 64000 -n 10 -R 'rusage[mem=64000]' "${mapper}/gsnap -A sam -B 5 -t 10 -Q -D ${gsnap_index} -s Mus_musculus_ERCC -d Mus_musculus_ERCC --gunzip $< --gunzip $(word 2, $^) | ${samtools}/samtools view -bS -o $@ -"

#### Read counting
.PHONY: read-count
read-count: ${read_count}


${results}/%.txt: ${cur_results}/%.bam
	${bsub} -M 8000 -R 'rusage[mem=8000]' "${samtools}/samtools view $< | ${htseq_count}/htseq-count -s no - ${merge}/Mus_musculus_ERCC.gtf > $@"

