### crossmap the benchmark mapping coordinates
chain=liftover.chain
CrossMap bam ${chain} MS5_${rep}.bam MapGlod_CP_${rep}.bam

### mapping accuracy
QUERY=${SAMPLE_ID}.sorted.rg.bam
python Mapping.accuracy.py -g MapGlod_CP_${rep}.bam.bam -q ${QUERY} -o graph_acc -t 10


### genotyping
REF=Sscrofa11.1_chr5.fa
gatk HaplotypeCaller --native-pair-hmm-threads 50 -R ${REF} -I ${SAMPLE_ID}.sorted.rg.bam -O ${SAMPLE_ID}.vcf.gz

### genotyping evaluation 
# benchmark
gatk SelectVariants -select-type SNP -V MS5.vcf.gz -O benchmark.SNP.vcf.gz
gatk SelectVariants -select-type INDEL -V MS5.vcf.gz -O benchmark.INDEL.vcf.gz
# evaluation
rtg format Sscrofa11.1_chr5.fa -o Scrofa11_chr5
VCF=Sus-pan_MC.vcf.gz
gatk SelectVariants -select-type SNP -V ${VCF} -O Pan.SNP.vcf.gz &
gatk SelectVariants -select-type INDEL -V ${VCF} -O Pan.INDEL.vcf.gz &
bench_SNP=benchmark.SNP.vcf.gz
bench_INDEL=benchmark.INDEL.vcf.gz
rtg vcfeval -b ${bench_SNP} -c Pan.SNP.vcf.gz -o Pan.SNP -t Scrofa11_chr5
rtg vcfeval -b ${bench_INDEL} -c Pan.INDEL.vcf.gz -o Pan.INDEL -t Scrofa11_chr5



