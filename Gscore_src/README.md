## Usage
### 1DEG_pearson_correlation.py
```shell
python 1DEG_pearson_correlation.py -e ./sample_data_input/example_GEM.txt -g ./sample_data_input/example_allDEG.txt -o ./sample_data_output/example_pcc.txt
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-e | input_GEM | The gene expression profile (.txt). | ./sample_data_input/example_GEM.txt
-g | input_DEG_list | The DEG list file, the DEG IDs are seperated by "\n". | ./sample_data_input/example_allDEG.txt
-o | output_path | The output file name for Pearson correlation network. | ./sample_data_output/example_pcc.txt

### 2Gscore.py
```shell
python 2Gscore.py -n ./sample_data_output/example_pcc.txt -g ./sample_data_input/example_allDEG.txt -s ./sample_data_input/example_GeneSet_KEGG_v102.txt -q ./sample_data_input/example_query_DEG.txt -t 0.7 -io ./sample_data_output/example_individual_DEG_result.txt -Lo ./sample_data_output/example_DEG_list_result.txt
```
Argument | Variable | Description | Default value
------------ | ------------- | ------------- | -------------
-n | input_network | The Pearson correlation network. | ./sample_data_output/example_pcc.txt
-g | input_DEG_list | The DEG list file, the DEG IDs are seperated by "\n". | ./sample_data_input/example_allDEG.txt
-s | input_geneset | The gene set collection file. | ./sample_data_input/example_GeneSet_KEGG_v102.txt
-q | input_query | The query DEG file, the DEG IDs are seperated by "\n". | ./sample_data_input/example_query_DEG.txt
-t | input_cutoff | The threshold for Pearson correlation coefficient. | 0.7
-io | output_path_single | The output file name for individual DEG result. | ./sample_data_output/example_individual_DEG_result.txt
-Lo | output_path_com | The output file name for DEG list result. | ./sample_data_output/example_DEG_list_result.txt

More detail information is available on the tutorial of the website: https://gscore.ibsb.nycu.edu.tw/tutorial.html

## Required Dependencies

* Python 3.x
* scipy.stats
* statsmodels.stats
* numpy
* pandas
