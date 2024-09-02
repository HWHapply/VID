# Viral Infection Detector: Ensemble Machine Learning for Predicting Viral Infection in Single-cell Transcriptomics of Virus-Induced Cancers


![ ](https://img.shields.io/badge/python-3.10-blue) ![ ](https://img.shields.io/badge/license-MIT-green) ![ ](https://img.shields.io/badge/R-4.4.1-blue) ![ ](https://img.shields.io/badge/conda-24.7.1-green)  

## Citation
Please cite: Wenhao Han, Jiahui Hu, Kane Toh Hui Chen. Viral Infection Detector: Ensemble Machine Learning for Predicting Viral Infection in Single-cell Transcriptomics of Virus-Induced Cancers, dd Month yyyy, PROTOCOL (Version 1) available at Protocol Exchange [paper_DOI]

## Table of Contents
* [Installation](#installation)
  * [Hardware](#hardware)
  * [Software](#software)
  * [Environment Setup](#environment-setup)
* [User Tutorial](#user-tutorial)
  * [Usage](#usage)
  * [Demo](#demo)
  * [Docker Run](#docker-run)

* [Expert Usage](#expert-usage)
  * [16S rRNA gene sequencing data of OSCC patients](#human-microbiome)
  * [Whole metagenomics data of Tara Ocean](#ocean-microbiome)
  * [Transcriptomics data of NASH patients](#human-transcriptome)


## Installation


### Hardware
This model is designed to be executed on standard computational hardware. While the process is feasible on commonly available systems, utilizing more robust computational resources can significantly enhance execution speed and efficiency.

Development and testing were carried out on different machines:

| Machine | CPU | Memory | OS | 
|--------|-----|--------------|---|
| Macbook Pro | 8-core Apple M1 chip | 16G | macOS Sonoma Version 14.6.1 |
| Linux WorkStation | 16-core Intel Xeon E5-2620 v4 braodwell-EP * 2  | 96G | Ubuntu 22.04 | 
| ThinkPad | ... | 16G | Windows 10 |


### Software
- R v4.4.1 or newer (https://www.r-project.org)
- Python3 v3.10.14 (https://www.python.org)
- Conda v14.7.1 (https://github.com/conda)

#### R packages
- Seurat (https://satijalab.org/seurat/)
- SeuratDisk (https://github.com/mojaveazure/seurat-disk)

#### Python packages
- pandas (https://pandas.pydata.org)
- NumPy (https://numpy.org/)
- Scipy (https://scipy.org/install/)
- scikit-learn (https://scikit-learn.org)
- boruta (https://github.com/scikit-learn-contrib/boruta_py)
- xgboost (https://xgboost.readthedocs.io/en/stable/install.html)
- pytorch (https://pytorch.org/get-started/locally/)
- skorch (https://github.com/skorch-dev/skorch)
- harmonypy (https://github.com/slowkow/harmonypy)
- anndata (https://anndata.readthedocs.io/en/latest/)
- tqdm (https://github.com/tqdm/tqdm)
- Matplotlib (https://matplotlib.org/)
- umap (https://umap-learn.readthedocs.io/en/latest/)
- seaborn (https://seaborn.pydata.org/)
  
#### Docker image
Above software list provides the minimal requirements for the complete execution of xMarkerFinder locally. Alternatively, we provide a ready-to-use Docker image, enabling users to skip the software installation and environment setup (https://hub.docker.com/r/tjcadd2022/xmarkerfinder). Additionally, an interactive JupyterHub server (https://mybinder.org/v2/gh/tjcadd2020/xMarkerFinder/HEAD) is also available.

### Environment setup
#### 1. Conda installation (version 24.7.1 is recommended)
Please install conda according to your OS from (https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

#### 2. Compiling from source code
Clone the source github repository to local machine:
```
git clone https://github.com/HWHapply/VID.git
```
Change to the work directory:
```
cd VID
```

#### 3. Create a conda environment called 'vid_env' for the VID running:
Choose the configuration file (vid_env_xxx.yml) to create the environment based on your OS:
```
conda create -f vid_env_{your_OS}.yml -n vid_env
```
Activate the environment:
```
conda activate vid_env
```
###### Remove the conda environment with commands below: 
Dectivate the VID environment:
```
conda deactivate vid_env
```
Remove the VID environment:
```
conda env remove -n vid_env
```

#### 4. Make the VID globally accessible:
Make the scripts executable:
```
chmod +x ./*
```
Create the solf link 'run_vid' of source script:
```
sudo ln -s $(realpath ./run_vid.sh) /usr/local/bin/run_vid
```
Execute with 'sudo' command if the 'Permission deny' error occurred. 

#### 5. Validate the environment setup:
Run the command below under any directory:
```
run_vid --help
```
You will get the help message below if the environment setup is successful:
```
usage: run_vid.py [-h] [--h5ad_dir H5AD_DIR] [--data_dir DATA_DIR] [--meta_dir META_DIR] [--output_dir OUTPUT_DIR] [--marker_dir MARKER_DIR]
                  [--feature_dir FEATURE_DIR] [--clinical_column CLINICAL_COLUMN] [--batch_column BATCH_COLUMN] [--sample_column SAMPLE_COLUMN]
                  [--test_ratio TEST_RATIO] [--num_split NUM_SPLIT] [--metamodel METAMODEL] [--threshold THRESHOLD] [--average AVERAGE]
                  [--random_state RANDOM_STATE] [--n_jobs N_JOBS] [--verbose VERBOSE]
```


#### Docker image setup
To simplify the setup process, we provide a Docker image that eliminates the need for manual equipment configuration. Please follow the steps below to use the Docker image:
- Install Docker:
Download and install Docker based on your operating system by following the instructions provided here: [Docker Installation Guide](https://docs.docker.com/engine/install/).
- Pull the Docker Image:
Retrieve the VID Docker image from Docker Hub using the following command:
```
docker pull hwhapply/vid:latest
```
- Verify the Image:
Confirm that the image has been successfully pulled by listing the available Docker images:
```
docker images
```
You should see 'hwhapply/vid:latest' listed under the REPOSITORY column.
- Execute Scripts within the Docker Container:
All subsequent scripts and procedures should be executed within a Docker container created from the hwhapply/vid image. The tutorial is in the next section.


## User tutorial
### Usage ###
Simply run VID in terminal with command below:
```
run_vid seuratobj_dir/xxx.rds
```
Please ensure that the input file conforms to the standard input format below. 
### __Input files__:  <br>
There are two input files that are required for VID running: <br>
#### 1. Seurat object (rds) ####
The input file seurat object saved in 'xxx.rds' format which generated with seurat single cell pipeline, the 'seuratobj_dir' is parent directory. 
```
seuratobj_dir/xxx.rds
```
Some columns should be included in metadata of the seurat object, visualize the metadata with code below:
```
# Run the code below in Rstudio or Jupyter:
library(seurat)
seurat_obj <- ReadRDS('seuratobj_dir/xxx.rds')
seurat_obj@@meta.data
```
The ideal metadata looks like the table below:
| orig.ident | ... | clinical_column |
|--------|-----|--------------|
| cell1_uid | ... | positive | 
| cell2_uid | ... |negative| 
| ... | ... | ... | ... |
| celln_uid |  ... | negative | 

- clinical_column ('clinical_column' by default): The sample level infection diagnosis, only has two str values: 'positive' and 'negative'.
- sample_column ('orig.ident' by default): orig.ident: The unique identifier for the cell origin (sample), included in the metadata by default.

You can also specify the corresponding names of those columns in your dataset with parameters 'clinical_column' and 'sample_column':
```
run_vid seuratobj_dir/xxx.rds \
        --clinical_column your_clinical_colname \
        --sample_column your_sample_id_colname
```
Specify the 'batch_column' as None if no batch appeared in dataset. The VID will ignore batch correction step if None is passed.

#### 2. Virus markers (txt) ####
A txt file contains the list of virus biomarkers, should be included in current work directory:
```
./markers.txt
```
Alternatively, specify your own virus marker file directory wtih parameter 'marker_dir':
```
run_vid seuratobj_dir/xxx.rds \
        --clinical_column your_clinical_colname \
        --sample_column your_sample_id_colname \
        --marker_dir your_marker_file_directory
```
### __Output files__: <br>
The code will automatically create output directory in current work directory named with the starting timestamp :
```
YYYYmmdd_HHMMSS 
├── data
│   ├── data.rds
│   └── metadata.csv
└── output
    ├── Confusion_Matrix_test.png
    ├── ROC_Curve_test.png
    ├── important_genes.txt
    ├── pred_proba_hist_test.png
    ├── pred_proba_hist_unseen.png
    ├── test_scores_weighted.csv
    ├── val_cv_scores_weighted.csv
    └── vid_YYmmdd_HHMMSS.pkl
```
#### Explanation: ####

- **YYYYmmdd_HHMMSS/**: The root directory, named with the current timestamp:

  - **data/**: Contains the input data and metadata table with results saved inside.
    - ***data.rds***: The input Seurat object file with predicted infection status and probabilities in the meta.data.
    - ***metadata.csv***: A CSV file containing metadata the same as the samples in `data.rds`.

  - **output/**: Contains the results and outputs from the machine learning tasks.
    - ***Confusion_Matrix_test.png***: An image file showing the confusion matrix on the test set.
    - ***ROC_Curve_test.png***: An image file showing the Receiver Operating Characteristic (ROC) curve on the test set.
    - ***important_genes.txt***: A text file listing the important genes identified by the boruta.
    - ***pred_proba_hist_test.png***: A histogram showing the distribution of predicted probabilities on the test dataset.
    - ***pred_proba_hist_unseen.png***: A histogram showing the distribution of predicted probabilities on an unseen dataset.
    - ***test_scores_weighted.csv***: A table containing the weighted test scores.
    - ***val_cv_scores_weighted.csv***: A table file containing the weighted cross-validation scores.
    - ***vid_YYmmdd_HHMMSS.pkl***: The VID object (for expert usage), timestamped with the current date and time.

### Parameters ###
__seuratobj_dir__ : str, requied
   > The directory of the input rds file (seurat object).

__output_dir__ : str, optional, default = './'
   > The output directory, set as current work directory by default.

__marker_dir__ : str, optional, default = ./markers.txt
   > The directory of a txt file contains the list of virus biomarkers, with each gene occupying one line.
   > The markers will be applied for the defination of traning set(truely infected and truely uninfected), 
   > while the cells from infected sample with any marker expressed will be considered truely infected.
   > The markers will also be ***excluded*** from modeling.

__feature_dir__ : str, optional, default = None
   > The directory of a txt file contains a list of important genes, with each gene occupying one line.
   > If given, ignore feature selection and apply the important genes for modeling, otherwise, the boruta
   > feature selection will be applied to select important genes.
  
__clinical_column__ : str, optional, default = clinical_column
   > The column indicates the sample-level clinical diagnosis of infection, the column should
   > included in the metadata of the input seurat object. 'positive' and 'negative' are two valid values in
   > this column. This column will be applied for training set defination together with 'marker_dir'.
   > Please prepare the column with valid values in metadata and pass the column name when running the code.

__batch_column__ : str, optional, default = batch
   > The column indicates the batch label which will be applied in batch correction with harmony.

__sample_column__ : str, optional, default = orig.ident
   > The column indicates the sample id. The column will be applied for training set defination.
   > Only the cell from uninfected sample without any virus markers expressed are considered to be truely uninfected.

__test_ratio__ : float, optional, default = 0.3
   > The ratio of test set. The testing set is splitted from the training data (infection status confirmed cells) before
   > model training. The test set will be applied as unseen data for model evaluation.

__num_split__ : int, optional, default = 5
   > The number of splitting for base model training and hyperparameter tuning of meta model.

__metamodel__ : str, optional {xgb, mlp}, default = xgb
   > The classifier applied as meta model. If 'xgb' passed, extreme gradient boosting classifier will be applied. If 'mlp'
   > passed, the multi-layer perceptron will be applied.

__threshold__ : float, optional {recommended range: 0.6 ~ 0.9}, default = None
   > The threshold for the decision function of final prediction. It can be understand as confidence of prediction: with higher threshold, 
   > the detected infection will be more reliable, but it may leads to misdetection of potential infected cells if the threshold is too high.
   > With this parameter specified, additional column with predicted infection status defined by this threshold will be added in meta data.
   > The default 'infection_status' in final meta data is defined with threshold = 0.5, the additional colunm defined by this threshold 
   > will be named as 'infection_status_{threshold}' which won't overwrite the the default. 

__average__ : str, optional {micro, macro, samples, weighted, binary, None} , default = weighted 
   > Define the type of averaging performed on the evaluation scores among different class. 

__random_state__ : int, optional , default = 42
   > The seed used to initialize a pseudo-random number generator (PRNG). It ensures that the results of random processes,
   > such as shuffling data or splitting datasets or model construction, are reproducible.

__n_jobs__ : int, optional, default = -1
   > The number of CPU applied for parallel excecution. All available CPU will be applied if not specify (n_jobs = -1).
   > Running with more CPUs can accelerate the training, but may effects other application and process on the computer.
   > Use this parameter with caution.

__verbose__ : int, optional, default = 2
   > Controls verbosity of output.

__help__ : Flag
   >  Show the help message and exit.

### Demo ###

#### NPC-EBV-Lymphocytes ####
Make a work directory for demo:
```
mkdir -p ./demo/data
```
Download the demo data to the demo directory from dropbox with wget:
```
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/bdkv2napos1md1uca2wg8/demo.rds?rlkey=bhe5deyz2o6kenj2s2fypxkzv&st=8armlfka&dl=1' -O ./demo/data/demo.rds
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/w4hhojmrnprcvj5emqljj/EBV_markers.txt?rlkey=rsa5lfha6rs8feivzmt5xpe9g&st=v13yukg0&dl=0' -O ./demo/data/EBV_markers.txt
```
You can also download and save data to './demo/data' with [demodata](https://www.dropbox.com/scl/fi/bdkv2napos1md1uca2wg8/demo.rds?rlkey=bhe5deyz2o6kenj2s2fypxkzv&st=8armlfka&dl=1).

Running VID:
```
run_vid ./demo/data/demo.rds \
        --output_dir ./demo \
        --marker_dir ./demo/data/EBV_markers.txt \
        --clinical_column ebv_status \
        --metamodel xgb 
```

#### NPC-EBV-Epithelial ####
Make a work directory for demo2:
```
mkdir -p ./demo2/data
```
Download the demo data to the demo directory from dropbox with wget:
```
wget --no-check-certificate 'https://www.dropbox.com/scl/fi/7lxap1hltlxqgy0v9vb05/demo2.rds?rlkey=eh4unrw8qkzz5wogl1ey3zozm&st=d91bko3t&dl=0' -O ./demo2/data/demo2.rds
```
You can also download and save data to './demo2/data' with [demodata2](https://www.dropbox.com/scl/fi/7lxap1hltlxqgy0v9vb05/demo2.rds?rlkey=eh4unrw8qkzz5wogl1ey3zozm&st=d91bko3t&dl=0).

Running VID:
```
run_vid ./demo2/data/demo2.rds \
        --output_dir ./demo2 \ 
        --marker_dir ./demo2/data/EBV_markers.txt \
        --clinical_column EBV_state \
        --metamodel mlp
```

### Docker Run ###

## Expert Usage ##
