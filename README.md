# taxdb-integration
format NCBI refseq genomes for non-prokaryotic taxa to integrate with GTDB

### Suggested workflow
* follow approach implemented in [kraken2](https://github.com/DerrickWood/kraken2) to access and download genomes for each non-prokaryotic section of refseq: fungi, invertebrate, plant, protozoa, vertebrate-mammalian, vertebrate-other, viral
* according to the NCBI readme: "The genomic.fna.gz file includes all top-level sequences in the assembly (chromosomes, plasmids, organelles, unlocalized scaffolds, unplaced scaffolds, and any alternate loci or patch scaffolds)", therefore the separate refseq plasmid, chloroplast, mitochondrion databases are not included additionally
* TODO: implement genome assessment for eukaryotic genomes from NCBI (e.g. [BUSCO](https://busco.ezlab.org/busco_userguide.html))
* dereplicate genomes per refseq section using approach suggested to [correct index databases](https://github.com/rrwick/Metagenomics-Index-Correction)
* for kraken2: parse taxonomy using [tax_from_gtdb.py](https://github.com/rrwick/Metagenomics-Index-Correction/blob/master/tax_from_gtdb.py) and directly build kraken2 database with all genomes (gtdb plus ncbi) in one go
* Alternative approach: merge with GTDB using [FlexTaxD](https://pypi.org/project/flextaxd/) **(not implemented)**
* TODO: format database for [kaiju](https://github.com/bioinformatics-centre/kaiju) (rely on available annotation for eukaryotic and viral genomes, update prokaryotic annotation)
