# AvrRmo1
Identification of *AvrRmo1*

## Motivation
To identify the gene underlying recognition by *Rmo1* in *Magnaporthe oryzae*. Mutagenesis and natural variation were used to identify *AvrRmo1*. UV mutagenesis was carried out on *M. oryzae* isolate Ken54-20. Approximately 15 independent mutants were identified.

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

data = read.table(file="Mo_Ken5420_jellyfish_24mer.histo.ID", header=T)
data = data.frame(data)

postscript(file="Mo_Ken5420_jellyfish_24mer_distribution.ps", width=6, height=4)
ggplot(data, aes(k, count)) + geom_point() + xlim(c(4,400)) + ylim(c(0,3e6)) + xlab("Frequency") + ylab("Total counts")
dev.off()
```

## Assembly of *M. oryzae* isolate Ken54-20

### Kmergenie

```bash
./kmergenie-1.7051/kmergenie DG_001_1.txt
```

> best k: 119

```bash
./kmergenie-1.7051/kmergenie DG_001_2.txt
```

> best k: 109


### Minia

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

The final assembly had the following statistics:

> N50 655,584
> Sequence 44,835,622
> Average 311,358
> E-size 748,645
> Count 144

### Assessment of assembly

#### Augustus

```bash
augustus --species=magnaporthe_grisea --strand=both GCF_000002495.2_MG8_genomic.fa > GCF_000002495.2_MG8_genomic.augustus.gff3 2>&1 &
augustus --species=magnaporthe_grisea --strand=both MoKen5420.masurca.contigs.fa > MoKen5420.masurca.contigs.augustus.gff3 2>&1 &
```

#### BUSCO

```bash
gffread GCF_000002495.2_MG8_genomic.augustus.gff3 -g GCF_000002495.2_MG8_genomic.fa -x GCF_000002495.2_MG8_genomic.augustus.cds.fa
gffread MoKen5420.masurca.contigs.augustus.gff3 -g MoKen5420.masurca.contigs.fa -x MoKen5420.masurca.contigs.augustus.cds.fa

python scripts/run_BUSCO.py -i GCF_000002495.2_MG8_genomic.augustus.cds.fa -l ascomycota_odb9 -o GCF_000002495.2_MG8_genomic.augustus.cds -m transcriptome -c 4
python scripts/run_BUSCO.py -i MoKen5420.masurca.contigs.augustus.cds.fa -l ascomycota_odb9 -o MoKen5420.masurca.contigs.augustus.cds -m transcriptome -c 4
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
> INFO	BUSCO analysis done. Total running time: 273.26234889 seconds
> INFO	Results written in /home/ubuntu/blast/busco/run_GCF_000002495.2_MG8_genomic.augustus.cds/

> M. oryzae Ken54-20
> INFO	Results:
> INFO	C:95.7%[S:95.2%,D:0.5%],F:2.1%,M:2.2%,n:1315
> INFO	1258 Complete BUSCOs (C)
> INFO	1252 Complete and single-copy BUSCOs (S)
> INFO	6 Complete and duplicated BUSCOs (D)
> INFO	28 Fragmented BUSCOs (F)
> INFO	29 Missing BUSCOs (M)
> INFO	1315 Total BUSCO groups searched



## Identifying SNPs in wild-type versus mutants

### BBmap

```bash
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.fa in1=DG_001_2_400_1.fq in2=DG_001_2_400_2.fq minid=0.98 maxindel=1 outm=MoKen5420_wt_mt_400.sam
bbmap/bbmap.sh ref=MoKen5420.masurca.contigs.fa in1=DG_001_2_600_1.fq in2=DG_001_2_600_2.fq minid=0.98 maxindel=1 outm=MoKen5420_wt_mt_600.sam
```

### VarScan

```bash
samtools view -@ 4 -f2 -Shub -o MoKen5420_wt_mt_400.bam MoKen5420_wt_mt_400.sam
samtools sort -@ 4 -o MoKen5420_wt_mt_400.sorted.bam MoKen5420_wt_mt_400.bam
samtools rmdup MoKen5420_wt_mt_400.sorted.bam MoKen5420_wt_mt_400.sorted.rmdup.bam

samtools view -@ 4 -f2 -Shub -o MoKen5420_wt_mt_600.bam MoKen5420_wt_mt_600.sam
samtools sort -@ 4 -o MoKen5420_wt_mt_600.sorted.bam MoKen5420_wt_mt_600.bam
samtools rmdup MoKen5420_wt_mt_600.sorted.bam MoKen5420_wt_mt_600.sorted.rmdup.bam

samtools merge MoKen5420_wt_mt.sorted.rmdup.bam MoKen5420_wt_mt_400.sorted.rmdup.bam MoKen5420_wt_mt_600.sorted.rmdup.bam

samtools index -@ 4 MoKen5420_wt_mt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_wt_mt.sorted.rmdup.bam > MoKen5420_wt_mt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_wt_mt.sorted.rmdup.pileup.txt > MoKen5420_wt_mt.sorted.rmdup.pileupindel.txt
```

```bash
samtools view -@ 4 -f2 -Shub -o MoKen5420_wt_wt_400.bam MoKen5420_wt_wt_400.sam
samtools sort -@ 4 -o MoKen5420_wt_wt_400.sorted.bam MoKen5420_wt_wt_400.bam
samtools rmdup MoKen5420_wt_wt_400.sorted.bam MoKen5420_wt_wt_400.sorted.rmdup.bam

samtools view -@ 4 -f2 -Shub -o MoKen5420_wt_wt_600.bam MoKen5420_wt_wt_600.sam
samtools sort -@ 4 -o MoKen5420_wt_wt_600.sorted.bam MoKen5420_wt_wt_600.bam
samtools rmdup MoKen5420_wt_wt_600.sorted.bam MoKen5420_wt_wt_600.sorted.rmdup.bam

samtools merge MoKen5420_wt_wt.sorted.rmdup.bam MoKen5420_wt_wt_400.sorted.rmdup.bam MoKen5420_wt_wt_600.sorted.rmdup.bam

samtools index -@ 4 MoKen5420_wt_wt.sorted.rmdup.bam
samtools mpileup -f MoKen5420.masurca.contigs.fa -BQ0 MoKen5420_wt_wt.sorted.rmdup.bam > MoKen5420_wt_wt.sorted.rmdup.pileup.txt
bedtools genomecov -d -split -ibam MoKen5420_wt_wt.sorted.rmdup.bam > MoKen5420_wt_wt.sorted.rmdup.genomecov.txt

java -jar VarScan.v2.3.8.jar mpileup2snp MoKen5420_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_wt_wt.sorted.rmdup.pileupsnp.txt
java -jar VarScan.v2.3.8.jar mpileup2indel MoKen5420_wt_wt.sorted.rmdup.pileup.txt > MoKen5420_wt_wt.sorted.rmdup.pileupindel.txt
```



### BBmap callvariants.sh

```bash
./bbmap/callvariants.sh in=MoKen5420_wt_mt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.fa out=MoKen5420_wt_mt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_wt_mt.sorted.rmdup.callvariants.gff3
./bbmap/callvariants.sh in=MoKen5420_wt_wt.sorted.rmdup.bam ref=MoKen5420.masurca.contigs.fa out=MoKen5420_wt_wt.sorted.rmdup.callvariants.vcf outgff=MoKen5420_wt_wt.sorted.rmdup.callvariants.gff3
```


```bash
gffread MoKen5420.masurca.contigs.augustus.gff3 -g MoKen5420.masurca.contigs.fa -y MoKen5420.masurca.contigs.augustus.protein.fa
python EffectorP.py -i ../../MoKen5420.masurca.contigs.augustus.protein.fa
```

```bash
python analyze_SNPs.py
```

```bash

```
