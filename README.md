# Effects of nitrogen starvation on gene expression in methylotrophic yeast Komagataella phaffii. Bioinformatic analysis of sequencing results

Authors:

- Makeeva Anastasiya
- Sidorin Anton
- Ishtuganova Valeria
- Sambuk Elena
- Padkina Marina
- Rumyantsev Andrey

Organization: St. Petersburg State University, 7/9 Universitetskaya nab., St. Petersburg, 199034 Russia.

## Libraries

Bash:

```
Trimmomatic=0.22
fastqc=0.11.9
gffread=0.12.1
hisat2=2.2.1
samtools=1.7
featureCounts=2.0.1
R=3.6.3
Rscript=3.6.3
python=3.8.12
```

R:

```
stringr (version: 1.4.0)
DESeq2 (version: 1.24.0)
EnhancedVolcano (version: 1.13.2)
dplyr (version: 1.0.8)
```

Python:

```
https://github.com/lbg-spbu/blast_Pp2Sc
```

## 1. Download raw data

Clone this repository and download raw data.

*TODO add link after publication in NCBI GEO*

```shell
git clone ...
cd 1.RawData
```

For further analysis the `1.RawData` folder should look like this:

```
1.RawData/
├── Control_1.fastq
├── Control_2.fastq
├── Control_3.fastq
├── Starv_1.fastq
├── Starv_2.fastq
└── Starv_3.fastq
```

## 2. Process raw data

We used **FastQC** ([link](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
program to evaluate the quality of reads and **Trimmomatic** ([doi](https://doi.org/10.1093/bioinformatics/btu170)) to
filter and trim reads.

```shell
cd ../2.TrimmedData
for i in {Control_1,Control_2,Control_3,Starv_1,Starv_2,Starv_3};
do echo "==================== $i ====================";
TrimmomaticSE -threads 4 ../1.RawData/"$i".fastq "$i".fastq \
  ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
  LEADING:25 \
  TRAILING:25 \
  SLIDINGWINDOW:4:25 \
  MINLEN:36 2>> stat.txt
echo "==================== end $i ====================";
done;

mkdir ./fastqc
fastqc *.fastq -o ./fastqc
```

## 3. Download reference genome data

To map trimmed reads to reference genome, we first need to obtain this
genome ([link](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000027005.1/)).

Also, since we need the GTF file for further analysis, we need to make it from GFF.

```shell
cd ../
mkdir P.p_genome && cd P.p_genome
curl -o P.p_genome.zip \
  'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000027005.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&hydrated=FULLY_HYDRATED'
unzip P.p_genome.zip
gffread ncbi_dataset/data/GCF_000027005.1/genomic.gff -T -o ncbi_dataset/data/GCF_000027005.1/genomic.gtf
```

## 4. Align trimmed reads to reference genome

We used **hisat-2** ([doi](https://doi.org/10.1038/s41587-019-0201-4)) to build genome indexes and perform alignment.

Also, we used **samtools** ([doi](https://doi.org/10.1093/gigascience/giab008)) for translation to BAM format.

```shell
mkdir index && cd index
hisat2-build -f \
  ../ncbi_dataset/data/GCF_000027005.1/GCF_000027005.1_ASM2700v1_genomic.fna \
  p_genome_index 2> Index_Logs.txt
```

```shell
cd ../../3.Alignment
for i in {Control_1,Control_2,Control_3,Starv_1,Starv_2,Starv_3};
do echo "==================== $i ====================";
hisat2 -x ../P.p_genome/index/p_genome_index -U ../2.TrimmedData/"$i".fastq -S "$i".sam 2>> stat.txt;
samtools view -b "$i".sam -o "$i".bam;
samtools sort "$i".bam > "$i".sorted.bam;
samtools index "$i".sorted.bam; 
rm "$i".bam "$i".sam;
echo "==================== end $i ====================";
done;
```

## 5. Count aligned reads

We used **featureCounts** ([doi](https://doi.org/10.1093/bioinformatics/btt656)) to count aligned reads.

```shell
cd ../4.TableCount
featureCounts \
 -g gene_id \
 -a ../P.p_genome/ncbi_dataset/data/GCF_000027005.1/genomic.gtf \
 -o Table_count.txt \
 ../3.Alignment/*.sorted.bam 2>> stat.txt
```

## 6. Differentially expressed gene analysis

We used **DESeq2** ([doi](https://doi.org/10.1186/s13059-014-0550-8)) to get differential gene expression.

```shell
cd ../5.DESeq
Rscript Pp_Nitrogen_Starv_DEG_Analysis.R
```

## 7. BLAST to _Saccharomyces cerevisiae_

To get information about differential expressed genes, we use our script, that makes a protein blast on _Saccharomyces
cerevisiae_. Thus, as a result, we get more information about each differentially expressed gene.

```shell
cd ../6.Blast_Pp2Sc

# 1. Download repo with code
git clone git@github.com:lbg-spbu/blast_Pp2Sc.git blast_script && cd blast_script
# 2. Create virtual environment
python3.8 -m venv venv
source venv/bin/activate
# 3. Download libraries from requirements.txt
pip install -r requirements.txt
# 4. Run script
python ./main.py -f ../../5.DESeq/Result_Starv_0.05_activated.tsv -o ../Pp2Sc_Starv_activated.tsv
python ./main.py -f ../../5.DESeq/Result_Starv_0.05_repressed.tsv -o ../Pp2Sc_Starv_repressed.tsv
# 5. Exit from virtual environment
deactivate

cd ..
rm -rf blast_script/
```
