## <a name="intro"></a>Introduction

Some organisms, such as paramecium, have MIC and MAC kernels. These kernels have different content and organisations. The MAC kernel in included in the larger MIC kernel. This tutorial is meant to show you how to assemble the MAC kernel for a dataset which is a mixture of MIC and MAC reads. 
First, we take as granted that MAC reads are many time more present in the read set than MIC reads (in some cases 100 times more present for example). Four time more present in the example used in this tutorial.  
Second MAC chromosomes which are many times more numerous as MIC chromosomes end with characteristic telomeric repeats (AAAACC and AAACCC or TTGGGG and TTTGGG in our example). These repeats will be indicators of correct placement of chromosome ends. 

## <a name="proc"></a>Procedure

The procedure uses kmer profiles to separate MIC and MAC reads. MAC reads should contain kmers with high coverage (arount the coverage mode of the complete read set) with only small low coverage sections corresponding to sequencing errors. MIC reads should contain much larger low coverage section corresponding to MIC specific regions which have not been amplified in the MIC to MAC transition. 

### <a name="readcov"></a>Kmer read coverage profile

The first step it to build the kmer read profile to find out at which coverage MIC and MAC kernels has been sequenced. This is performed using jellyfish to create a kmer dictionary and genomescope2 to generate the profile graphical view. 

We will use the public read set ERR11843474 which can be downloaded using the following link 
[ERR11843474]: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR118/074/ERR11843474/ERR11843474.fastq.gz

To install jellyfish and genomescope2 see : 
[jellyfish]: https://github.com/gmarcais/Jellyfish
[genomescope2]: https://github.com/schatzlab/genomescope

```sh
# read to kmer histogram file 
jellyfish count -m 21 -C -o hifi.jf -s 100M -t 24 <(gzip -d -c ERR11843474.fastq.gz)
jellyfish histo -h 1000000 hifi.jf > hifi.hist
Rscript genomescope.R -i hifi.hist -k 21 -p 2 -o genomescope2
```

Genomescope2 resulting image 
![Read set kmer profile](https://github.com/chklopp/macassemblies/blob/main/transformed_linear_plot.png)

We consider the peak at 190x coverage to be MIC specific kmer and the coverages above to be MAC related. 

Now we know that the MAC kernel coverage mode is 829x.

We can now generate a kmer coverage profile along each read. This profile will be used to select the MAC reads. 
These profiles are generated using the jellyfish query_per_sequence script which can be found here : 
[query_per_sequence]: https://github.com/gmarcais/Jellyfish/blob/master/examples/query_per_sequence/query_per_sequence.cc

But before we have to convert our read fastq  file in fasta

```sh
gunzip -c ERR11843474.fastq.gz | sed -n '1~4s/^@/>/p;2~4p' > ERR11843474.fasta
```

Now to generate the kmer profiles 
```sh
query_per_sequence hifi.jf ERR11843474.fasta > ERR11843474.fasta.kmer_profiles
```

### <a name="profilefilt"></a>Kmer profile filtering to keep only MAC reads

To filter MAC reads we have to find the size of low coverage sections in the reads. The kmer_profiles file contains the kmer coverage for each position in the reads (but the 20 last postions). We are going to calculate the median coverage for each 50 base pairs window using the median_split.awk script. Using a kmer of 21 the median over 50 bases pairs should not decrease for nucleotide errors in reads but only for MIC sections. 

```sh
cat ERR11843474.fasta.kmer_profiles | paste - - | sed 's/>//' | awk -f median_split.awk > ERR11843474.fasta.kmer_profiles.median
```

From this new file we can calculate the number or chunks having a median lower than 10 and the number of chunks having a median over 1000. This is performed using and awk command line.

```sh
awk '{a=0;b=0;for (i=2;i<=NF;i++){if ($i < 250){a++}; if ($i > 2500){b++}} print $1"\t"a"\t"b}' ERR11843474.fasta.kmer_profiles.median > ERR11843474.fasta.kmer_profiles.median.stats
```

We remove from the list all the reads having at least one section with a coverage under 10. 

```sh
awk '$2 == 0{print $1}'   ERR11843474.fasta.kmer_profiles.median.stats >  ERR11843474.fasta.kmer_profiles.median.stats.tokeep 
samtools faidx ERR11843474.fasta
samtools faidx ERR11843474.fasta ERR11843474.fasta.kmer_profiles.median.stats.tokeep > MAC_reads.fasta
```

### <a name="assembly"></a>MAC reads assembly with hifiasm
We then assemble MAC reads using hifiasm 

```sh
hifiasm -t 16 --hg-size 100m  -o hifiasm_MAC MAC_reads.fasta
awk '/^S/{print ">"$2"\n"$3}' hifiasm_MAC.bp.p_ctg.gfa | fold > hifiasm_MAC.bp.p_ctg.gfa.fa
awk '/^S/{print $2"\t"$4"\t"$5}' hifiasm_MAC.bp.p_ctg.gfa | sed 's/LN:i://;s/rd:i://' > hifiasm_MAC.bp.p_ctg.gfa.cov
```

### <a name="assemblyfilter"></a>MAC assembly filtering 

```sh
hifiasm -t 16 --hg-size 100m  -o hifiasm_MAC MAC_reads.fasta
```


### <a name="assemblycheck"></a>MAC assembly check 
