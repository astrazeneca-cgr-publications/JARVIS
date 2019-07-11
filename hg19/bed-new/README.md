# Downlod Human GRCh37 gtf from Ensembl (same as gencode hg19)
wget ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz


# Get merged Exon regions
zcat < Homo_sapiens.GRCh37.87.chr.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sortBed | mergeBed -i - | gzip > Homo_sapiens.GRCh37.87.Exons.merged.bed.gz

# Get merged Intron regions
zcat < Homo_sapiens.GRCh37.87.chr.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | subtractBed -a stdin -b Homo_sapiens.GRCh37.87.Exons.merged.bed.gz | gzip > Homo_sapiens.GRCh37.87.Introns.merged.bed.gz

# Get Intergenic regions
mysql --user=genome --host=genome-euro-mysql.soe.ucsc.edu -A -e "select chrom, size from hg19.chromInfo"  > hg19.genome

zcat < Homo_sapiens.GRCh37.87.chr.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg19.genome | gzip > Homo_sapiens.GRCh37.87.Intergenic.merged.bed.gz
