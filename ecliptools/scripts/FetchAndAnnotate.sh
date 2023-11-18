# Get HepG2 Cell Line DHX30 eCLIP IDR Data
wget https://www.encodeproject.org/files/ENCFF663QIZ/@@download/ENCFF663QIZ.bed.gz
gunzip ENCFF663QIZ.bed.gz

# Get K562 Cell Line DHX30 eCLIP IDR Data
wget https://www.encodeproject.org/files/ENCFF128AKC/@@download/ENCFF128AKC.bed.gz
gunzip ENCFF128AKC.bed.gz

# Sort by chromosome then position using BEDTOOLS bed-sort
sort-bed ENCFF128AKC.bed > ENCFF128AKC.sorted.bed
sort-bed ENCFF128AKC.bed > ENCFF128AKC.sorted.bed

# Get ucsc whole genome comprehensive annotation data
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/wgEncodeGencodeCompV43.txt.gz
gunzip wgEncodeGencodeCompV43.txt.gz

# Convert into .bed file
java -jar ../jvarkit/dist/kg2bed.jar wgEncodeGencodeCompV43.txt > knownGenes.bed
sort-bed knownGenes.bed > knownGenes.sorted.bed

# For later clarity append a column marking all rRNA locis with "rRNA".
# Also rearranges the columns to fit the formatting of knownGenes.bed
awk 'BEGIN {FS="\t"; OFS="\t"} {print $0, "rRNA"}' rRNA_loci.bed |\
 awk 'BEGIN {FS="\t"; OFS="\t"} {print $1, $2, $3, $6, $4, $7, $5}' > temp.bed
mv rRNA_loci.bed original.rRNA_loci.bed
mv temp.bed rRNA_loci.bed
sort-bed rRNA_loci.bed > rRNA_loci.sorted.bed

# Merge the knownGenes and rRNA annotation bed files for easier processing
# NOTE: those bed files are NOT .bed compliant they also contain some data delimited by ; instead of \t
cat knownGenes.sorted.bed rRNA_loci.sorted.bed > genesAndrRNA.bed
sort-bed genesAndrRNA.bed > genesAndrRNA.sorted.bed

# annotate the file using the BEDTOOLS bedmap utility
bedmap --delim "\t" --multidelim "\t" --unmapped-val "UNKNOWN" --echo --echo-map ENCFF663QIZ.sorted.bed genesAndrRNA.sorted.bed > ENCFF663QIZ.mapped.bed
bedmap --delim "\t" --multidelim "\t" --unmapped-val "UNKNOWN" --echo --echo-map ENCFF128AKC.sorted.bed genesAndrRNA.sorted.bed > ENCFF128AKC.mapped.bed


