## Transcriptomics-based Network Medicine Discovery and Population-based Validation Identifies Metformin as a Candidate Drug for Atrial Fibrillation

### Network proximity code
The proximity code can be found in the `proximity/` folder. First, decompress the file `HumanInteractome.7z` in the same folder.

* HumanInteractome.tsv - the human interactome, including sources and evidence types
* HumanInteractome.npy - numpy matrix of all precalculated shortest distance in the interactome
* DrugTargetNetwork.txt - drug target network
* network_proximity.py - code to compute closest network proximity. See below

The `network_proximity.py` supports two modes. To run the program (Python 3), numpy and networkx need to be installed.
```
pip install numpy networkx
```

#### To compute the closest proximity between two gene lists
```
python network_proximity.py path_to_gene_list_1 path_to_gene_list_2 number_of_repeat random_seed
```

#### To screen all the drugs for a gene list
```
python network_proximity.py DRUG path_to_gene_list number_of_repeat random_seed
```

#### Gene list file format
Save the gene Entrez IDs in a single column, e.g.,
```
58
162
166
310
373
396
409
427
476
481
514
```
