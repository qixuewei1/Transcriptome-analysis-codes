# QC
ngsqc -1 1_S179_L002_R1_001.fastq.gz -2 1_S179_L002_R2_001.fastq.gz -o 01.fastq-qc -k HZ:1
Rscript ngsqc.R --base HZ:1.atgc --qual HZ:1.qual --out 01.fastq-qc/fig/HZ:1
ngsqc -1 2_S180_L002_R1_001.fastq.gz -2 2_S180_L002_R2_001.fastq.gz -o 01.fastq-qc -k JA:1
Rscript ngsqc.R --base JA:1.atgc --qual JA:1.qual --out 01.fastq-qc/fig/JA:1
Rscript $Bin/bin/ngsqc.R --base list.atgc --qual outlist.qual --out sample.result.list

# assemble
cat 1_S179_L002_R1_001.fastq.gz 2_S180_L002_R1_001.fastq.gz > 02.assemble/read1.fq.gz
cat 1_S179_L002_R2_001.fastq.gz 2_S180_L002_R2_001.fastq.gz > 02.assemble/read2.fq.gz
Trinity --seqType fq --left read1.fq.gz --right read2.fq.gz --max_memory 200G --min_contig_length 200 --CPU 12 --min_kmer_cov 2 --max_chrysalis_cluster_size 20 --min_per_id_same_path 98  --max_memory 150G --output Trinity_outdir
perl lenth.pl -i Trinity.fasta -o 02.assemble/Trinity.fasta -len 200
hisat2-build -p 8 02.assemble/Trinity.fasta 02.assemble/Trinity.fasta

# transrate
hisat2 -p 8 --dta  -x 02.assemble/Trinity.fasta -1 read1.fq.gz -2 read2.fq.gz -S mapping.sam  2 >mapping.stat
transrate --assembly=02.assemble/Trinity.fasta --left=read1.fq.gz --right=read2.fq.gz --threads=8 --output=03.transrate2
cd outdir/ && TransDecoder.LongOrfs -t ref.fasta --gene_trans_map mapping.stat --output_dir outdir/TransDecoder
TransDecoder.Predict -t ref.fasta -T 3000 --output_dir outdir/TransDecoder
Trinity_gene_splice_modeler.py  --trinity outdir/ref.trans.fa --out_prefix outdir/ref.genome

# predict
ln -s outdir/ref.trans.fa.transdecoder.gff3 outdir/ref.trans.gff
gffread ref.gff -g ref.fa -T -o outdir/ref.trans.gtf -x outdir/ref.trans.cds.fa  -y outdir/ref.trans.pro.fa
extract-transcript-to-gene-map-from-trinity outdir/ref.trans.fa outdir/knowIso.form.txt
rsem-prepare-reference --bowtie2 -p 8 --gtf outdir/ref.trans.gtf --transcript-to-gene-map map.result outdir/ref.trans.fa outdir/ref.trans.fa
extract_exons.py outdir/ref.genome.gtf > outdir/ref.genome.exon && extract_splice_sites.py outdir/ref.genome.gtf > outdir/ref.genome.ss
hisat2-build -p 8 outdir/ref.genome.fasta --ss outdir/ref.genome.ss --exon outdir/ref.genome.exon outdir/ref.genome.fasta
samtools faidx outdir/ref.genome.fasta
samtools dict outdir/ref.genome.fasta -o outdir/ref.genome.dict

#express
rsem-calculate-expression -p 8 --bowtie2 --paired-end  fq1.list fq2.list ref.fasta  outdir/
cd $outdir/ && rsem-generate-data-matrix genes.list >outdir/genes.express.matrix && sed -i 's/\\.genes\\.results//g' $outdir/genes.express.matrix && sed -i 's/\"//g' outdir/genes.express.matrix
cd $outdir/ && rsem-generate-data-matrix isoforms.list> outdir/isoforms.express.matrix && sed -i 's/\\.isoforms\\.results//g' outdir/isoforms.express.matrix &&  sed -i 's/\"//g' outdir/isoforms.express.matrix
cd $outdir/ && Rscript density_ggplot2.R -m outdir/genes.express.matrix -o outdir && Rscript density_ggplot2.R -m outdir/isoforms.express.matrix -o outdir
cd $outdir/ && Rscript express_venn.R  -m outdir/genes.express.matrix -o outdir -d genes && Rscript express_venn.R -m outdir/isoforms.express.matrix -o outdir -d isoforms

#diff
Rscript matrix.R --infile matrix --outfile outdir/ --group grouplist
python split_group_contrix.py --grouplist  grouplist --compair  compare  --outdir outdir
Rscript DESeq2_v4.R -m genes.list -g outdir/*._group.list -o outdir/genes 
Rscript DESeq2_v4.R -m isoforms.list -g outdir/*.group.list -o outdir/isoforms 
python arrange_dexp22delist.py -i outdir/gene_diff.list -o outdir -f genes_dill.list
python arrange_dexp22delist.py -i outdir/isoforms_diff.list -o outdir -f isoforms_dill.list
Rscript difflist_venn.R -m $outdir/genes_dill.list -o outdir -d genes
Rscript difflist_venn.R -m $outdir/isoforms_dill.list -o outdir -d isoforms 

