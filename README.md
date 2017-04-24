# _WAS
Runs both Genome Wide Association Study (GWAS) and Exome Wide Association Study (EWAS) on VCF format genotype data using 24 different combinations of frequency, function, and relative filters. Used in O'Beirne et al (in review).

Software Dependencies:
PLINK https://www.cog-genomics.org/plink2
SnpEff http://snpeff.sourceforge.net/
VCFTools http://vcftools.sourceforge.net/
CADD http://cadd.gs.washington.edu/
R https://cran.r-project.org/
SKAT https://cran.r-project.org/web/packages/SKAT/index.html
EMMAX http://genetics.cs.ucla.edu/emmax/
KING http://people.virginia.edu/~wc9c/KING/
Java http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
Python https://www.python.org/download/releases/3.0/
Bgzip/Tabix http://www.htslib.org/download/ 


Data Dependencies:
Input genotype data in VCF format https://samtools.github.io/hts-specs/VCFv4.2.pdf
Input phenotype and covariate data in PLINK format https://www.cog-genomics.org/plink/1.9/input#covar
1000 Genomes site list with allele frequency in VCF format. ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
ExAC site list with allele frequency in VCF format. ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz

