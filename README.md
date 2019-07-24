# AvrRmo1
Identification of *AvrRmo1*

## Motivation
To identify the gene underlying recognition by *Rmo1* in *Magnaporthe oryzae*. Mutagenesis and natural variation were used to identify *AvrRmo1*. UV mutagenesis was carried out on *M. oryzae* isolate Ken54-20. Approximately 15 independent mutants were identified.

## Trimmomatic

```bash
java -jar trimmomatic-0.39.jar PE -threads 16 -phred33 DG_001_1_400_1.fq DG_001_1_400_2.fq DG_001_1_400_gDNA_forward_paired.fq DG_001_1_400_gDNA_forward_unpaired.fq DG_001_1_400_gDNA_reverse_paired.fq DG_001_1_400_gDNA_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > DG_001_1_400_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.39.jar PE -threads 16 -phred33 DG_001_1_600_1.fq DG_001_1_600_2.fq DG_001_1_600_gDNA_forward_paired.fq DG_001_1_600_gDNA_forward_unpaired.fq DG_001_1_600_gDNA_reverse_paired.fq DG_001_1_600_gDNA_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > DG_001_1_600_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.39.jar PE -threads 16 -phred33 DG_001_2_400_1.fq DG_001_2_400_2.fq DG_001_2_400_gDNA_forward_paired.fq DG_001_2_400_gDNA_forward_unpaired.fq DG_001_2_400_gDNA_reverse_paired.fq DG_001_2_400_gDNA_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > DG_001_2_400_trimmomatic.run.log 2>&1 &
java -jar trimmomatic-0.39.jar PE -threads 16 -phred33 DG_001_2_600_1.fq DG_001_2_600_2.fq DG_001_2_600_gDNA_forward_paired.fq DG_001_2_600_gDNA_forward_unpaired.fq DG_001_2_600_gDNA_reverse_paired.fq DG_001_2_600_gDNA_reverse_unpaired.fq ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:10 MINLEN:36 > DG_001_2_600_trimmomatic.run.log 2>&1 &
```


## Jellyfish
To confirm the quality of Illumina sequencing and confirm genome size of *M. oryzae* isolate Ken54-20, we used `jellyfish`.

```bash
jellyfish count -t 72 -C -m 17 -s 30G -o MoKen5420_jellyfish_17mer DG_001_1_400_1.fq DG_001_1_400_2.fq DG_001_1_600_1.fq DG_001_1_600_2.fq
jellyfish histo -h 3000000 -o MoKen5420_jellyfish_17mer.histo MoKen5420_jellyfish_17mer

jellyfish count -t 72 -C -m 24 -s 30G -o MoKen5420_jellyfish_24mer DG_001_1_400_1.fq DG_001_1_400_2.fq DG_001_1_600_1.fq DG_001_1_600_2.fq
jellyfish histo -h 3000000 -o MoKen5420_jellyfish_24mer.histo MoKen5420_jellyfish_24mer

jellyfish count -t 48 -C -m 17 -s 30G -o MoKen5420m1_jellyfish_17mer DG_001_2_400_1.fq DG_001_2_400_2.fq DG_001_2_600_1.fq DG_001_2_600_2.fq
jellyfish histo -h 3000000 -o MoKen5420m1_jellyfish_17mer.histo MoKen5420m1_jellyfish_17mer

jellyfish count -t 48 -C -m 24 -s 30G -o MoKen5420m1_jellyfish_24mer DG_001_2_400_1.fq DG_001_2_400_2.fq DG_001_2_600_1.fq DG_001_2_600_2.fq
jellyfish histo -h 3000000 -o MoKen5420m1_jellyfish_24mer.histo MoKen5420m1_jellyfish_24mer
```

Both wild-type and mutant had identical *k*-mer distributions that fit the expected model of a haploid species.

> GenomeScope version 1.0
> k = 24
> 
> property                      min               max               
> Heterozygosity                0.0595689%        0.0607236%        
> Genome Haploid Length         49,919,238 bp     49,927,451 bp     
> Genome Repeat Length          12,870,673 bp     12,872,791 bp     
> Genome Unique Length          37,048,565 bp     37,054,661 bp     
> Model Fit                     95.4049%          96.6713%          
> Read Error Rate               0.199458%         0.199458%         

To visualise the results in `R`, the following code can be used after adding a header to the `jellyfish histo` output.

```R
library(ggplot2)

data = read.table(file="MoKen5420_jellyfish_17mer.histo.ID", header=T)
data = data.frame(data)

postscript(file="MoKen5420_jellyfish_17mer_distribution.ps", width=6, height=4)
ggplot(data, aes(k, count)) + geom_point() + xlim(c(6,400)) + ylim(c(0,1e6)) + xlab("Frequency") + ylab("Total counts")
dev.off()

data = read.table(file="MoKen5420m1_jellyfish_17mer.histo.ID", header=T)
data = data.frame(data)

postscript(file="MoKen5420m1_jellyfish_17mer_distribution.ps", width=6, height=4)
ggplot(data, aes(k, count)) + geom_point() + xlim(c(6,400)) + ylim(c(0,1e6)) + xlab("Frequency") + ylab("Total counts")
dev.off()
```

## Oxford Nanopore Technologies sequencing
Oxford Nanopore Technologies (ONT) sequencing was used on *M. oryzae* isolate Ken54-20. A single flowcell was used for sequencing based on two runs. Guppy (v3.2.2) was used for base calling of raw data from ONT sequencing.

```bash
guppy_basecaller -i mgg_ken54 -s . --kit SQK-LSK109 --flowcell FLO-MIN106 -r --compress_fastq --cpu_threads_per_caller 96 --num_callers 1 > guppy_basecaller.log 2>&1 &
guppy_basecaller -i mgg_ken54_20b -s mgg_ken54_20b_results --kit SQK-LSK109 --flowcell FLO-MIN106 -r --compress_fastq --cpu_threads_per_caller 96 --num_callers 1 > guppy_basecaller_2.log 2>&1 &
```

Read quality assessment was performed using [pauvre](https://github.com/conchoecia/pauvre) and high quality reads identified and trimmed using [NanoFilt](https://github.com/wdecoster/nanofilt/).

```bash
pauvre stats --fastq mgg_ken54_20.fastq > mgg_ken54_20.fastq.pauvre.stats 2>&1 &
pauvre stats --fastq mgg_ken54_20b.fastq > mgg_ken54_20b.fastq.pauvre.stats 2>&1 &

cat mgg_ken54_20.fastq mgg_ken54_20b.fastq > mgg_ken54_20duo.fastq
pauvre stats --fastq mgg_ken54_20duo.fastq > mgg_ken54_20duo.fastq.pauvre.stats 2>&1 &
gzip -c9 mgg_ken54_20duo.fastq

gunzip -c fastq_runid_72ad5b4694bd5d06ce9c9be9a574a74b8f384a20_0.fastq.gz | NanoFilt -q 12 --headcrop 75 | gzip > trimmed_0b.fastq.gz.fastq.gz
```

NanoPore sequencing result for June 2019.

> numReads: 848,839
> %totalNumReads: 100.00
> numBasepairs: 11,257,175,152
> %totalBasepairs: 100.00
> meanLen: 13,262
> medianLen: 4,529
> minLen: 1
> maxLen: 238,793
> N50: 35,770
> L50: 94,871

## Assembly of *M. oryzae* isolate Ken54-20
### Kmergenie
`Kmergenie` was used to evaluate potential *k*-mers to use for `minia` assembly.

```bash
./kmergenie-1.7051/kmergenie DG_001_1.txt
```

> best k: 119

```bash
./kmergenie-1.7051/kmergenie DG_001_2.txt
```

> best k: 109


### Minia
*De novo* assembly of the genome of *M. oryzae* isolate Ken54-20 and mutant 1 is first performed using `minia`. This will provide a resource for future evaluation of the hybrid assembly and *AvrRmo1* candidate gene analysis. In this case, we used *k* = 119 for both assemblies for consistency.

```bash
./minia-v3.2.1-bin-Linux/bin/minia -in DG_001_1.txt -nb-cores 72 -kmer-size 119 -out MoKen5420.minia.k119
./minia-v3.2.1-bin-Linux/bin/minia -in DG_001_2.txt -nb-cores 72 -kmer-size 119 -out MoKen5420m1.minia.k119
```

### MaSuRCA
The hybrid assembler [MaSuRCA](https://github.com/alekseyzimin/masurca/) was used to assemble the genome of *M. oryzae* isolate Ken54-20. Version 3.3.3 was used for assembly. Default parameters were used throughout. The `sr_config.txt` file was modified to include Nanopore and Illumina data for this project. The commands are relatively straight forward. 

```bash
./install.sh
./assembly.sh
```

We built three assemblies (v1, v2, and v3) using an original poor quality ONT sequencing reaction (v1), raw reads from the recent ONT sequencing run (v2), and trimmed, quality controlled reads from the recent ONT sequencing run (v3). The latter assembly (v3) was identified as the best due to a range of parameters that will be discussed at length below.

The final assembly had the following statistics:

> v1
> N50 655,584
> Sequence 44,835,622
> Average 311,358
> E-size 748,645
> Count 144

> v2
> N50 5,941,966
> Sequence 47,427,932
> Average 1.05395e+06
> E-size 7.17211e+06
> Count 45

> v3
> N50 3,032,866
> Sequence 47,904,193
> Average 1.22831e+06
> E-size 3.62593e+06
> Count 39

#### Genome polishing
[Pilon](https://github.com/broadinstitute/pilon) was used to correct the *M. oryzae* Ken54-20 genome. Alignment of reads was performed using Bowtie2 and BWA, with BWA selected as the preferred aligner due to lower called SNPs and manually assessment of putative SNPs (i.e. BWA has less false positives).

##### Bowtie2
```bash
bowtie2-build --threads 16 MoKen5420.masurca.contigs.v3.fa MoKen5420.masurca.contigs.v3

bowtie2 -p 16 -x MoKen5420.masurca.contigs.v3 -1 DG_001_1_400_1.fq,DG_001_1_600_1.fq -2 DG_001_1_400_2.fq,DG_001_1_600_2.fq -S MoKen5420_wt_wt_v3.sam

samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_wt_v3.bam MoKen5420_wt_wt_v3.sam
samtools sort -@ 16 -o MoKen5420_wt_wt_v3.sorted.bam MoKen5420_wt_wt_v3.bam
samtools rmdup MoKen5420_wt_wt_v3.sorted.bam MoKen5420_wt_wt_v3.sorted.rmdup.bam

./bbmap/callvariants.sh in=MoKen5420_wt_wt_v3.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.v3.fa vcf=MoKen5420_wt_wt_v3.sorted.rmdup.vcf

samtools index -@ 16 MoKen5420_wt_wt_v3.sorted.rmdup.bam

java -Xmx100G -jar pilon-1.23.jar --genome MoKen5420.masurca.contigs.v3.fa --frags MoKen5420_wt_wt_v3.sorted.rmdup.bam --output MoKen5420_v3_pilon
```

##### BWA-MEM
```bash
bwa index MoKen5420.masurca.contigs.v3.fa

bwa mem -t 96 MoKen5420.masurca.contigs.v3.fa DG_001_1_400_1.fq DG_001_1_400_2.fq > MoKen5420_wt_wt_v3_400.sam
bwa mem -t 96 MoKen5420.masurca.contigs.v3.fa DG_001_1_600_1.fq DG_001_1_600_2.fq > MoKen5420_wt_wt_v3_600.sam

samtools view -@ 96 -f2 -Shub -o MoKen5420_wt_wt_v3_400.bam MoKen5420_wt_wt_v3_400.sam
samtools sort -@ 96 -o MoKen5420_wt_wt_v3_400.sorted.bam MoKen5420_wt_wt_v3_400.bam
samtools rmdup MoKen5420_wt_wt_v3_400.sorted.bam MoKen5420_wt_wt_v3_400.sorted.rmdup.bam

samtools view -@ 96 -f2 -Shub -o MoKen5420_wt_wt_v3_600.bam MoKen5420_wt_wt_v3_600.sam
samtools sort -@ 96 -o MoKen5420_wt_wt_v3_600.sorted.bam MoKen5420_wt_wt_v3_600.bam
samtools rmdup MoKen5420_wt_wt_v3_600.sorted.bam MoKen5420_wt_wt_v3_600.sorted.rmdup.bam

samtools merge -@ 96 MoKen5420_wt_wt_v3.sorted.rmdup.bam MoKen5420_wt_wt_v3_400.sorted.rmdup.bam MoKen5420_wt_wt_v3_600.sorted.rmdup.bam

./bbmap/callvariants.sh in=MoKen5420_wt_wt_v3.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.v3.fa vcf=MoKen5420_wt_wt_v3.sorted.rmdup.vcf

samtools index -@ 96 MoKen5420_wt_wt_v3.sorted.rmdup.bam

java -Xmx100G -jar pilon-1.23.jar --genome MoKen5420.masurca.contigs.v3.fa --frags MoKen5420_wt_wt_v3.sorted.rmdup.bam --output MoKen5420_v3_pilon
```

### Assessment of assembly

#### K-mer Analysis Toolkit (KAT)
We assessed the quality 

[KAT](https://github.com/TGAC/KAT)

```bash
kat hist -t 96 DG_001_1_*
kat gcp -t 96 DG_001_1_*
kat comp -t 96 'DG_001_1_400_1.fq DG_001_1_400_2.fq DG_001_1_600_1.fq DG_001_1_600_2.fq' MoKen5420.masurca.contigs.v3.4.fa
```

#### Augustus
*Ab initio* gene prediction was carried out using Augustus (v3.3.2) with the *M. oryzae* species gene model prediction.

```bash
augustus --species=magnaporthe_grisea --strand=both GCF_000002495.2_MG8_genomic.fa > GCF_000002495.2_MG8_genomic.augustus.gff3 2>&1 &
augustus --species=magnaporthe_grisea --strand=both MoKen5420.masurca.contigs.fa > MoKen5420.masurca.contigs.augustus.gff3 2>&1 &
augustus --species=magnaporthe_grisea --strand=both MoKen5420.masurca.contigs.v2.fa > MoKen5420.masurca.contigs.v2.augustus.gff3 2>&1 &
augustus --species=magnaporthe_grisea --strand=both MoKen5420.masurca.contigs.v2.pilon.fa > MoKen5420.masurca.contigs.v2.pilon.augustus.gff3 2>&1 &
augustus --species=magnaporthe_grisea --strand=both MoKen5420.masurca.contigs.v3.4.fa > MoKen5420.masurca.contigs.v3.4.augustus.gff3 2>&1 &
```

#### BUSCO

```bash
gffread GCF_000002495.2_MG8_genomic.augustus.gff3 -g GCF_000002495.2_MG8_genomic.fa -x GCF_000002495.2_MG8_genomic.augustus.cds.fa
gffread MoKen5420.masurca.contigs.augustus.gff3 -g MoKen5420.masurca.contigs.fa -x MoKen5420.masurca.contigs.augustus.cds.fa
gffread MoKen5420.masurca.contigs.v2.augustus.gff3 -g MoKen5420.masurca.contigs.v2.fa -x MoKen5420.masurca.contigs.v2.augustus.cds.fa
gffread MoKen5420.masurca.contigs.v2.pilon.augustus.gff3 -g MoKen5420.masurca.contigs.v2.pilon.fa -x MoKen5420.masurca.contigs.v2.pilon.augustus.cds.fa
gffread MoKen5420.masurca.contigs.v3.4.augustus.gff3 -g MoKen5420.masurca.contigs.v3.4.fa -x MoKen5420.masurca.contigs.v3.4.augustus.cds.fa

python scripts/run_BUSCO.py -i GCF_000002495.2_MG8_genomic.augustus.cds.fa -l ascomycota_odb9 -o GCF_000002495.2_MG8_genomic.augustus.cds -m transcriptome -c 4
python scripts/run_BUSCO.py -i MoKen5420.masurca.contigs.augustus.cds.fa -l ascomycota_odb9 -o MoKen5420.masurca.contigs.augustus.cds -m transcriptome -c 4
python scripts/run_BUSCO.py -i MoKen5420.masurca.contigs.v2.augustus.cds.fa -l ascomycota_odb9 -o MoKen5420.masurca.contigs.v2.augustus.cds -m transcriptome -c 16
python scripts/run_BUSCO.py -i ../annotation/MoKen5420.masurca.contigs.v2.pilon.augustus.cds.fa -l ascomycota_odb9 -o MoKen5420.masurca.contigs.v2.pilon.augustus.cds -m transcriptome -c 16
python scripts/run_BUSCO.py -i ../annotation/MoKen5420.masurca.contigs.v3.4.augustus.cds.fa -l ascomycota_odb9 -o MoKen5420.masurca.contigs.v3.4.augustus.cds -m transcriptome -c 16
```

> M. oryzae 70-15
> INFO	Results:
> INFO	C:97.2%[S:97.0%,D:0.2%],F:2.3%,M:0.5%,n:1315
> INFO	1279 Complete BUSCOs (C)
> INFO	1276 Complete and single-copy BUSCOs (S)
> INFO	3 Complete and duplicated BUSCOs (D)
> INFO	30 Fragmented BUSCOs (F)
> INFO	6 Missing BUSCOs (M)
> INFO	1315 Total BUSCO groups searched

> M. oryzae Ken54-20 v1
> INFO	Results:
> INFO	C:95.7%[S:95.2%,D:0.5%],F:2.1%,M:2.2%,n:1315
> INFO	1258 Complete BUSCOs (C)
> INFO	1252 Complete and single-copy BUSCOs (S)
> INFO	6 Complete and duplicated BUSCOs (D)
> INFO	28 Fragmented BUSCOs (F)
> INFO	29 Missing BUSCOs (M)
> INFO	1315 Total BUSCO groups searched

> M. oryzae Ken54-20 v2
> INFO	Results:
> INFO	C:98.0%[S:96.1%,D:1.9%],F:1.7%,M:0.3%,n:1315
> INFO	1289 Complete BUSCOs (C)
> INFO	1264 Complete and single-copy BUSCOs (S)
> INFO	25 Complete and duplicated BUSCOs (D)
> INFO	22 Fragmented BUSCOs (F)
> INFO	4 Missing BUSCOs (M)
> INFO	1315 Total BUSCO groups searched

> M. oryzae Ken54-20 v3.4
> INFO	Results:
> INFO  C:97.9%[S:95.5%,D:2.4%],F:1.7%,M:0.4%,n:1315
> INFO  1288    Complete BUSCOs (C)
> INFO  1256    Complete and single-copy BUSCOs (S)
> INFO  32      Complete and duplicated BUSCOs (D)
> INFO  23      Fragmented BUSCOs (F)
> INFO  4       Missing BUSCOs (M)
> INFO  1315    Total BUSCO groups searched

## Evaluating the assembly using long NanoPore reads

```bash
gunzip -c mgg_ken54_20duo.fastq.gz | NanoFilt -q 12 --headcrop 75 -l 10000 | gzip > mgg_ken54_20duo_trimmed.fastq.gz
minimap2 -ax map-ont MoKen5420.masurca.contigs.v2.fa ../nanopore/mgg_ken54_20duo_trimmed.fastq.gz > MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.sam
samtools view -@ 16 -F4 -Shub -o MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.bam MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.sam
samtools sort -@ 16 -o MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.sorted.bam MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.bam
bedtools genomecov -d -split -ibam MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.sort.bam > MoKen5420.masurca.contigs.v2_mgg_ken54_20duo_trimmed.sort.genomecov.txt

minimap2 -ax map-ont MoKen5420.masurca.contigs.v2.pilon.fa nanopore/mgg_ken54_20duo_trimmed.fastq.gz > MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.sam
samtools view -@ 16 -F4 -Shub -o MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.bam MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.sam
samtools sort -@ 16 -o MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.sorted.bam MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.bam
bedtools genomecov -d -split -ibam MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.sorted.bam > MoKen5420.masurca.contigs.v2.pilon_mgg_ken54_20duo_trimmed.sorted.genomecov.txt

minimap2 -ax map-ont MoKen5420.masurca.contigs.v3.4.fa mgg_ken54_20duo_trimmed.fastq -t 72 > MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.sam
samtools view -@ 72 -F4 -Shub -o MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.bam MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.sam
samtools sort -@ 72 -o MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.sorted.bam MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.bam
bedtools genomecov -d -split -ibam MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.sorted.bam > MoKen5420.masurca.contigs.v3.4_mgg_ken54_20duo_trimmed.sorted.genomecov.txt
```

## Identifying SNPs in wild-type versus mutants

### BWA-MEM
```bash
bwa index MoKen5420.masurca.contigs.v3.4.fa

bwa mem -t 72 MoKen5420.masurca.contigs.v3.4.fa DG_001_1_400_1.fq DG_001_1_400_2.fq > MoKen5420_wt_wt_v3.4_400.sam
bwa mem -t 72 MoKen5420.masurca.contigs.v3.4.fa DG_001_1_600_1.fq DG_001_1_600_2.fq > MoKen5420_wt_wt_v3.4_600.sam

samtools view -@ 72 -f2 -Shub -o MoKen5420_wt_wt_v3.4_400.bam MoKen5420_wt_wt_v3.4_400.sam
samtools sort -@ 72 -o MoKen5420_wt_wt_v3.4_400.sorted.bam MoKen5420_wt_wt_v3.4_400.bam
samtools rmdup MoKen5420_wt_wt_v3.4_400.sorted.bam MoKen5420_wt_wt_v3.4_400.sorted.rmdup.bam

samtools view -@ 72 -f2 -Shub -o MoKen5420_wt_wt_v3.4_600.bam MoKen5420_wt_wt_v3.4_600.sam
samtools sort -@ 72 -o MoKen5420_wt_wt_v3.4_600.sorted.bam MoKen5420_wt_wt_v3.4_600.bam
samtools rmdup MoKen5420_wt_wt_v3.4_600.sorted.bam MoKen5420_wt_wt_v3.4_600.sorted.rmdup.bam

samtools merge -@ 72 MoKen5420_wt_wt_v3.4.sorted.rmdup.bam MoKen5420_wt_wt_v3.4_400.sorted.rmdup.bam MoKen5420_wt_wt_v3.4_600.sorted.rmdup.bam

bwa mem -t 72 MoKen5420.masurca.contigs.v3.4.fa DG_001_2_400_1.fq DG_001_2_400_2.fq > MoKen5420_wt_mt_v3.4_400.sam
bwa mem -t 72 MoKen5420.masurca.contigs.v3.4.fa DG_001_2_600_1.fq DG_001_2_600_2.fq > MoKen5420_wt_mt_v3.4_600.sam

samtools view -@ 72 -f2 -Shub -o MoKen5420_wt_mt_v3.4_400.bam MoKen5420_wt_mt_v3.4_400.sam
samtools sort -@ 72 -o MoKen5420_wt_mt_v3.4_400.sorted.bam MoKen5420_wt_mt_v3.4_400.bam
samtools rmdup MoKen5420_wt_mt_v3.4_400.sorted.bam MoKen5420_wt_mt_v3.4_400.sorted.rmdup.bam

samtools view -@ 72 -f2 -Shub -o MoKen5420_wt_mt_v3.4_600.bam MoKen5420_wt_mt_v3.4_600.sam
samtools sort -@ 72 -o MoKen5420_wt_mt_v3.4_600.sorted.bam MoKen5420_wt_mt_v3.4_600.bam
samtools rmdup MoKen5420_wt_mt_v3.4_600.sorted.bam MoKen5420_wt_mt_v3.4_600.sorted.rmdup.bam

samtools merge -@ 72 MoKen5420_wt_mt_v3.4.sorted.rmdup.bam MoKen5420_wt_mt_v3.4_400.sorted.rmdup.bam MoKen5420_wt_mt_v3.4_600.sorted.rmdup.bam
```

## Code below is depreciated
### BBmap

```bash
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.fa in1=DG_001_2_400_1.fq in2=DG_001_2_400_2.fq minid=0.98 maxindel=1 outm=MoKen5420_wt_mt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.fa in1=DG_001_2_600_1.fq in2=DG_001_2_600_2.fq minid=0.98 maxindel=1 outm=MoKen5420_wt_mt_600.sam

bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_1_400_1.fq in2=DG_001_1_400_2.fq minid=0.98 maxindel=1 outm=MoKen5420_v2_wt_wt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_1_600_1.fq in2=DG_001_1_600_2.fq minid=0.98 maxindel=1 outm=MoKen5420_v2_wt_wt_600.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_2_400_1.fq in2=DG_001_2_400_2.fq minid=0.98 maxindel=1 outm=MoKen5420_v2_wt_mt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_2_600_1.fq in2=DG_001_2_600_2.fq minid=0.98 maxindel=1 outm=MoKen5420_v2_wt_mt_600.sam

bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_1_400_1.fq in2=DG_001_1_400_2.fq outm=MoKen5420_v2_wt_wt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_1_600_1.fq in2=DG_001_1_600_2.fq outm=MoKen5420_v2_wt_wt_600.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_2_400_1.fq in2=DG_001_2_400_2.fq outm=MoKen5420_v2_wt_mt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.v2.fa in1=DG_001_2_600_1.fq in2=DG_001_2_600_2.fq outm=MoKen5420_v2_wt_mt_600.sam
```

### Bowtie2

```bash
bowtie2-build --threads 16 MoKen5420.masurca.contigs.v2.fa MoKen5420.masurca.contigs

bowtie2 -p 16 -x MoKen5420.masurca.contigs -1 DG_001_1_400_1.fq,DG_001_1_600_1.fq -2 DG_001_1_400_2.fq,DG_001_1_600_2.fq -S MoKen5420_wt_wt.sam
bowtie2 -p 16 -x MoKen5420.masurca.contigs -1 DG_001_2_400_1.fq,DG_001_2_600_1.fq -2 DG_001_2_400_2.fq,DG_001_2_600_2.fq -S MoKen5420_wt_mt.sam

samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_wt.bam MoKen5420_wt_wt.sam
samtools sort -@ 16 -o MoKen5420_wt_wt.sorted.bam MoKen5420_wt_wt.bam
samtools rmdup MoKen5420_wt_wt.sorted.bam MoKen5420_wt_wt.sorted.rmdup.bam

samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_mt.bam MoKen5420_wt_mt.sam
samtools sort -@ 16 -o MoKen5420_wt_mt.sorted.bam MoKen5420_wt_mt.bam
samtools rmdup MoKen5420_wt_mt.sorted.bam MoKen5420_wt_mt.sorted.rmdup.bam


samtools index -@ 16 MoKen5420_wt_mt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupindel.txt

```

### VarScan

```bash
samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_mt_400.bam MoKen5420_wt_mt_400.sam
samtools sort -@ 16 -o MoKen5420_wt_mt_400.sorted.bam MoKen5420_wt_mt_400.bam
samtools rmdup MoKen5420_wt_mt_400.sorted.bam MoKen5420_wt_mt_400.sorted.rmdup.bam

samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_mt_600.bam MoKen5420_wt_mt_600.sam
samtools sort -@ 16 -o MoKen5420_wt_mt_600.sorted.bam MoKen5420_wt_mt_600.bam
samtools rmdup MoKen5420_wt_mt_600.sorted.bam MoKen5420_wt_mt_600.sorted.rmdup.bam

samtools merge MoKen5420_wt_mt.sorted.rmdup.bam MoKen5420_wt_mt_400.sorted.rmdup.bam MoKen5420_wt_mt_600.sorted.rmdup.bam

samtools index -@ 16 MoKen5420_wt_mt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupindel.txt
```

```bash
samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_wt_400.bam MoKen5420_wt_wt_400.sam
samtools sort -@ 16 -o MoKen5420_wt_wt_400.sorted.bam MoKen5420_wt_wt_400.bam
samtools rmdup MoKen5420_wt_wt_400.sorted.bam MoKen5420_wt_wt_400.sorted.rmdup.bam

samtools view -@ 16 -f2 -Shub -o MoKen5420_wt_wt_600.bam MoKen5420_wt_wt_600.sam
samtools sort -@ 16 -o MoKen5420_wt_wt_600.sorted.bam MoKen5420_wt_wt_600.bam
samtools rmdup MoKen5420_wt_wt_600.sorted.bam MoKen5420_wt_wt_600.sorted.rmdup.bam

samtools merge MoKen5420_wt_wt.sorted.rmdup.bam MoKen5420_wt_wt_400.sorted.rmdup.bam MoKen5420_wt_wt_600.sorted.rmdup.bam

samtools index -@ 16 MoKen5420_wt_wt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_wt_wt.sorted.rmdup.bam > MoKen5420_wt_wt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_wt_wt.sorted.rmdup.bam > MoKen5420_wt_wt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_wt_wt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_wt_wt.sorted.rmdup.pileupindel.txt
```

```bash
samtools view -@ 16 -f2 -Shub -o MoKen5420_v2_wt_mt_400.bam MoKen5420_v2_wt_mt_400.sam
samtools sort -@ 16 -o MoKen5420_v2_wt_mt_400.sorted.bam MoKen5420_v2_wt_mt_400.bam
samtools rmdup MoKen5420_v2_wt_mt_400.sorted.bam MoKen5420_v2_wt_mt_400.sorted.rmdup.bam

samtools view -@ 16 -f2 -Shub -o MoKen5420_v2_wt_mt_600.bam MoKen5420_v2_wt_mt_600.sam
samtools sort -@ 16 -o MoKen5420_v2_wt_mt_600.sorted.bam MoKen5420_v2_wt_mt_600.bam
samtools rmdup MoKen5420_v2_wt_mt_600.sorted.bam MoKen5420_v2_wt_mt_600.sorted.rmdup.bam

samtools merge MoKen5420_v2_wt_mt.sorted.rmdup.bam MoKen5420_v2_wt_mt_400.sorted.rmdup.bam MoKen5420_v2_wt_mt_600.sorted.rmdup.bam

samtools index -@ 16 MoKen5420_v2_wt_mt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_v2_wt_mt.sorted.rmdup.bam > MoKen5420_v2_wt_mt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_v2_wt_mt.sorted.rmdup.bam > MoKen5420_v2_wt_mt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_v2_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_v2_wt_mt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_v2_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_v2_wt_mt.sorted.rmdup.pileupindel.txt
```

```bash
samtools view -@ 16 -f2 -Shub -o MoKen5420_v2_wt_wt_400.bam MoKen5420_v2_wt_wt_400.sam
samtools sort -@ 16 -o MoKen5420_v2_wt_wt_400.sorted.bam MoKen5420_v2_wt_wt_400.bam
samtools rmdup MoKen5420_v2_wt_wt_400.sorted.bam MoKen5420_v2_wt_wt_400.sorted.rmdup.bam

samtools view -@ 16 -f2 -Shub -o MoKen5420_v2_wt_wt_600.bam MoKen5420_v2_wt_wt_600.sam
samtools sort -@ 16 -o MoKen5420_v2_wt_wt_600.sorted.bam MoKen5420_v2_wt_wt_600.bam
samtools rmdup MoKen5420_v2_wt_wt_600.sorted.bam MoKen5420_v2_wt_wt_600.sorted.rmdup.bam

samtools merge MoKen5420_v2_wt_wt.sorted.rmdup.bam MoKen5420_v2_wt_wt_400.sorted.rmdup.bam MoKen5420_v2_wt_wt_600.sorted.rmdup.bam

samtools index -@ 16 MoKen5420_v2_wt_wt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_v2_wt_wt.sorted.rmdup.bam > MoKen5420_v2_wt_wt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_v2_wt_wt.sorted.rmdup.bam > MoKen5420_v2_wt_wt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_v2_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_v2_wt_wt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_v2_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_v2_wt_wt.sorted.rmdup.pileupindel.txt
```



### BBmap callvariants.sh

```bash
./bbmap/callvariants.sh in=MoKen5420_wt_mt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.fa out=MoKen5420_wt_mt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_wt_mt.sorted.rmdup.callvariants.gff3
./bbmap/callvariants.sh in=MoKen5420_wt_wt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.fa out=MoKen5420_wt_wt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_wt_wt.sorted.rmdup.callvariants.gff3

./bbmap/callvariants.sh in=MoKen5420_v2_wt_mt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.v2.fa out=MoKen5420_v2_wt_mt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_v2_wt_mt.sorted.rmdup.callvariants.gff3
./bbmap/callvariants.sh in=MoKen5420_v2_wt_wt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.v2.fa out=MoKen5420_v2_wt_wt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_v2_wt_wt.sorted.rmdup.callvariants.gff3
```


```bash
gffread MoKen5420.masurca.contigs.augustus.gff3 -g MoKen5420.masurca.contigs.fa -y MoKen5420.masurca.contigs.augustus.protein.fa
python EffectorP.py -i ../../MoKen5420.masurca.contigs.augustus.protein.fa

gffread MoKen5420.masurca.contigs.v2.augustus.gff3 -g MoKen5420.masurca.contigs.v2.fa -y MoKen5420.masurca.contigs.v2.augustus.protein.fa
python EffectorP.py -i ../../MoKen5420.masurca.contigs.v2.augustus.protein.fa

gffread MoKen5420.masurca.contigs.v2.pilon.augustus.gff3 -g MoKen5420.masurca.contigs.v2.pilon.fa -y MoKen5420.masurca.contigs.v2.pilon.augustus.protein.fa
python EffectorP.py -i ../../MoKen5420.masurca.contigs.v2.pilon.augustus.protein.fa
```

```bash
python analyze_SNPs.py
```

```bash
nucmer --maxmatch -l 100 -c 500 Maggr1_AssemblyScaffolds.fasta MoKen5420.masurca.contigs.v3.4.fa --prefix Maggr1_Ken5420 -t 4
nucmer --maxmatch -l 100 -c 500 Magnaporthe_oryzae.MG8.dna.toplevel.fa MoKen5420.masurca.contigs.v3.4.fa --prefix MG8_Ken5420 -t 4
nucmer --maxmatch -l 100 -c 500 guy11_smartdenovo_polished.fasta MoKen5420.masurca.contigs.v3.4.fa --prefix Guy11_Ken5420 -t 4
nucmer --maxmatch -l 100 -c 500 KE002_Pacbio_42_contigs.fa MoKen5420.masurca.contigs.v3.4.fa --prefix KE002_Ken5420 -t 4
```
