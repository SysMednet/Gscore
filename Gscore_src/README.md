## Usage
### 1DEG_pearson_correlation.py
```shell
python 1DEG_pearson_correlation.py -e ./sample_data_input/example_GEM.txt -g ./sample_data_input/example_allDEG.txt -o ./sample_data_output/example_pcc.txt
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-e | input_GEM | Input file of the gene expression data (.txt). | ./sample_data_input/example_GEM.txt
-g | input_DEG_list | Input file of the DEG list, in which Entrez IDs are separated by "\n". | ./sample_data_input/example_allDEG.txt
-o | output_path | Output file of the dataset-derived coexpression network. | ./sample_data_output/example_pcc.txt

### 2Gscore.py
```shell
python 2Gscore.py -n ./sample_data_output/example_pcc.txt -g ./sample_data_input/example_allDEG.txt -s ./sample_data_input/example_GeneSet_KEGG_v102.txt -q ./sample_data_input/example_query_DEG.txt -t 0.7 -io ./sample_data_output/example_individual_DEG_result.txt -Lo ./sample_data_output/example_DEG_list_result.txt
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-n | input_network | Input file of the dataset-derived coexpression network. | ./sample_data_output/example_pcc.txt
-g | input_DEG_list | Input file of the DEG list, in which Entrez IDs are separated by "\n". | ./sample_data_input/example_allDEG.txt
-s | input_geneset | Input file the collection of gene sets. | ./sample_data_input/example_GeneSet_KEGG_v102.txt
-q | input_query | Input file of the query DEG list, in which Entrez IDs are separated by "\n". | ./sample_data_input/example_query_DEG.txt
-t | input_cutoff | Criterion of |Pearson correlation coefficient| for determining co-expressed DEG pairs between the query  list and gene sets. | 0.7
-io | output_path_single | Output file of the analysis results of individual DEG . | ./sample_data_output/example_individual_DEG_result.txt
-Lo | output_path_com | Output file of the analysis results of the DEG list. | ./sample_data_output/example_DEG_list_result.txt

More detailed information is available at the Gscore website tutorial: https://gscore.ibsb.nycu.edu.tw/tutorial.html

## Required Dependencies

* Python 3.x
* SciPy
* statsmodels
* numpy
* pandas
