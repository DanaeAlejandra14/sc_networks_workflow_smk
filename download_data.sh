# For downloading data, activate the synapse.yml env with conda env create --file=workflow/envs/synapse.yaml --name=synapseclient
# conda activate synapseclient
#The tree must look like this

#data

#├── data
#│   ├── 190403-B4-A.barcodes.tsv.gz
#│   ├── 190403-B4-A.features.tsv.gz
#│   ├── 190403-B4-A.matrix.mtx.gz
#│   ├── 190403-B4-B.barcodes.tsv.gz
#│   ├── 190403-B4-B.features.tsv.gz
#│   ├── 190403-B4-B.matrix.mtx.gz
#│   ├── 191122-B6-R1710143-alone.barcodes.tsv.gz
#│   ├── 191122-B6-R1710143-alone.features.tsv.gz
#│   ├── 191122-B6-R1710143-alone.matrix.mtx.gz
#...

#resources (annotation and metadata)

#├── resources
#│   ├── Adult_Human_PrefrontalCortex.pkl
#│   ├── ROSMAP_biospecimen_metadata.csv
#│   ├── ROSMAP_clinical.csv
#│   ├── ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv
#│   ├── cell-annotation.full-atlas.csv
#│   ├── individual_phenotype_labels.csv
#│   └── tools

#Download metadata into resources/

#clinical metadata
#synapse get syn3191087 --downloadLocation resources/

#biospecimen
#synapse get syn21323366 --downloadLocation resources/

#demultiplex data
#synapse get syn34572333 --downloadLocation resources/

#cell annotation
#synapse get syn55219673 --downloadLocation resources/

#count matrices
#synapse get -r syn51123521 --downloadLocation data/

#Cell type annotation from celltypist

#wget https://celltypist.cog.sanger.ac.uk/models/Human_PFC_Ma/v1/Adult_Human_PrefrontalCortex.pkl
