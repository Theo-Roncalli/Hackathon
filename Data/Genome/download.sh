#!/usr/bin/env bash

echo "Downloading genome..."
wget "http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
echo "Done"

echo "Downloading annotations..."
wget "http://labshare.cshl.edu/shares/gingeraslab/www-data/dobin/STAR/STARgenomes/Human/GRCh38_Ensembl99_sparseD3_sjdbOverhang99/Homo_sapiens.GRCh38.99.gtf"
echo "Done"
