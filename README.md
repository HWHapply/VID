# Viral Infection Detector: Ensemble Learning for Predicting Viral Infection in Single-cell Transcriptomics of Virus-Induced Cancers

![ ](https://img.shields.io/badge/python-3.10-blue) ![ ](https://img.shields.io/badge/license-MIT-green) ![ ](https://img.shields.io/badge/R-4.4.1-blue) ![ ](https://img.shields.io/badge/conda-24.7.1-green)  ![ ](https://img.shields.io/badge/docker-27.1.1-blue)  

## Citation
Please cite: Wenhao Han, Jiahui Hu, Kane Toh, Hui Chen*. Viral Infection Detector: Ensemble Learning for Predicting Viral Infection in Single-cell Transcriptomics of Virus-Induced Cancers, dd Month yyyy, PROTOCOL (Version 1) available at Protocol Exchange [paper_DOI]

## Table of Contents
* [Installation](#installation)
  * [Hardware](#hardware)
  * [Software](#software)
  * [Environment Setup](#environment-setup)
* [User Tutorial](#user-tutorial)
  * [Usage](#usage)
  * [Demo](#demo)
  * [Docker Run](#docker-run)
* [Transfer VID Model On Unseen Data](#transfer-vid-model-on-unseen-data)
  * [Transfer Learning With Conda](#transfer-learning-with-conda)
  * [Transfer Learning With Docker](#transfer-learning-with-docker)

## Installation


### Hardware
This model is designed to be executed on standard computational hardware. While the process is feasible on commonly available systems, utilizing more robust computational resources can significantly enhance execution speed and efficiency.

Development and testing were carried out on different machines:

| Machine | CPU | Memory | OS | 
|--------|-----|--------------|---|
| Macbook Pro | 8-core Apple M1 chip | 16G | macOS Sonoma Version 14.6.1 |
| Linux WorkStation | 16-core Intel Xeon E5-2620 v4 braodwell-EP * 2  | 96G | Ubuntu 22.04 | 
| Alienware 17 R5 | 6-core Intel i7-8750H CPU | 16G | Windows 10 |



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
- seaborn (https://seaborn.pydata.org/)
  
#### Docker image
Above software list displays the minimal requirements for running the VID locally. Alternatively, to enabla users to skip the software installation and environment setup, an Docker image is provided [hwhapply/vid:latest](https://hub.docker.com/r/hwhapply/vid). 

### Environment setup
#### 1. Conda installation (version 24.7.1 is recommended)
Please install conda according to your OS from (https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

If you are windows user, please install git bash for your system from(https://git-scm.com/downloads/win). After installation, Press `Win + X`, then search `git bash` to open the git bash command line. Then proceed to the next step.

#### 2. Compiling from source code
Clone the source github repository to local machine:
```
git clone https://github.com/HWHapply/VID.git
```
Change to the work directory:
```
cd VID
```
Make sure the `setup.sh` is excutable and run it from `VID` directory:
```bash
chmod +x setup.sh
./setup.sh
```
#### 3. Validate the environment setup:
Activate the `vid_env`:
```bash
conda activate vid_env
```
Run the command below under any directory:
```bash
run_vid --help
```
You will get the help message below if the environment setup is successful:
```
Usage: /Your/local/path/to/run_vid   
    seuratobj_dir, required
                        The directory of input seurat object.
    --marker_dir MARKER_DIR, -mkd MARKER_DIR, required
                        The markers stores in a txt file(one gene per row).
    --clinical_column CLINICAL_COLUMN, -cc CLINICAL_COLUMN, required
                        The column indicates the infection status in clinical assessment.(Sample level)
    --output_dir OUTPUT_DIR, -od OUTPUT_DIR, optional
                        The output directory.
    --feature_dir FEATURE_DIR, -fd FEATURE_DIR, optional
                        The directory of txt file stores the important features(gene).
    --label_dir LABEL_DIR, -ld LABEL_DIR, optional
                        The directory of txt file stores the pre-defined labels.
    --batch_column BATCH_COLUMN, -bc BATCH_COLUMN, optional
                        The column indicates the batch label that will be used for batch correction(harmony).
    --sample_column SAMPLE_COLUMN, -sc SAMPLE_COLUMN, optional
                        The column indicates the sample id.
    --test_ratio TEST_RATIO, -tr TEST_RATIO, optional
                        The ratio of validating set.
    --num_split NUM_SPLIT, -ns NUM_SPLIT, optional
                        The number of splitting for base model training and hyperparameter tuning for meta model.
    --metamodel METAMODEL, -mm METAMODEL, optional
                        The classifier applied as meta model.
    --threshold THRESHOLD, -threds THRESHOLD, optional
                        The threshold for the decision function of final prediction.
    --average AVERAGE, -avg AVERAGE, optional
                        Define the type of averaging performed on the evaluation scores among different class.
    --random_state RANDOM_STATE, -rs RANDOM_STATE, optional
                        The random state for the reproduction of result.
    --n_jobs N_JOBS, -threads N_JOBS, optional
                        Number of threads applied for parallel excecution.
    --verbose VERBOSE, -v VERBOSE, optional
                        The verbose mode.
    --vidmodel_dir VIDMODEL_DIR, -vmd VIDMODEL_DIR, optional
                        The directory of the vid object to applied for transfer learning.
    -h, --help            
                        Show this help message and exit.
```
Please be aware that it may take few minutes for vid to install all required packages in the first running.

#### Docker image setup
To simplify the setup process, we provide a Docker image that eliminates the requirement for manual environment configuration. Please follow the steps below to use the Docker image:
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
The 'hwhapply/vid:latest' repository has listed under the REPOSITORY column.
```
REPOSITORY           TAG       IMAGE ID       CREATED         SIZE
hwhapply/vid         latest    a59825e4b92d   21 hours ago    7.96GB
```
- Execute Scripts within the Docker Container:
All subsequent scripts and procedures should be executed within a Docker container created from the 'hwhapply/vid' image. The tutorial is in the next section.


## User tutorial
### Usage ###
Simply run VID in terminal with command below:
```
run_vid seuratobj_dir/xxx.rds --marker_dir markers.txt --clinical_column clinical_column
```
Please ensure that the input file conforms to the standard input format specified below. 
### __Input files__:  <br>
There are two input files that are required for VID running: <br>
#### 1. Seurat object (rds) ####
The input file seurat object saved in 'xxx.rds' format which generated with seurat single cell pipeline, the `seuratobj_dir` is parent directory. 
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

- `clinical_column` ('clinical_column' by default): The sample level infection status, only has two str values: 'positive' and 'negative'.
- `sample_column` ('orig.ident' by default): The unique identifier for sample the cell originate from, included in the metadata by default.

You can also specify those columns in your dataset accordingly with parameters `clinical_column` and `sample_column`:
```
run_vid seuratobj_dir/xxx.rds \
--clinical_column your_clinical_colname \
--sample_column your_sample_id_colname
```
Specify the `batch_column` as None if no batch effect appeared in dataset. The VID will ignore batch correction step if None is passed.

#### 2. Viral markers (txt) ####
A text file contains the list of viral biomarkers, should be specified with parameter `marker_dir`, the content of marker file shows below:
```
viral marker 1
viral marker 2
viral marker 3
...
```
The markers will be filetered out from the high variance genes before model training. If the label file is not specified with parameter `label_dir` by user, the markers will also be applied for labeling.

### __Output files__: <br>
The code will automatically create output directory in current working directory named with the starting timestamp:
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
    ├── cv_box.png
    ├── vid_YYmmdd_HHMMSS.pkl
    └── log.txt

```
#### Explanation: ####

- **YYYYmmdd_HHMMSS**: The root directory, named with the current timestamp:

  - **data**: Contains the input data and metadata table with results saved inside.
    - ***data.rds***: The input Seurat object file with predicted infection status and probabilities in the meta.data.
    - ***metadata.csv***: A CSV file containing metadata the same as the samples in `data.rds`.

  - **output**: Contains the results and outputs from the machine learning tasks.
    - ***Confusion_Matrix_test.png***: An image file showing the confusion matrix on the test set.
    - ***ROC_Curve_test.png***: An image file showing the Receiver Operating Characteristic (ROC) curve on the test set.
    - ***important_genes.txt***: A text file listing the important genes identified by the boruta.
    - ***pred_proba_hist_test.png***: A histogram showing the distribution of predicted probabilities on the test dataset.
    - ***pred_proba_hist_unseen.png***: A histogram showing the distribution of predicted probabilities on an unseen dataset.
    - ***test_scores_weighted.csv***: A table containing the weighted test scores.
    - ***val_cv_scores_weighted.csv***: A table file containing the weighted cross-validation scores.
    - ***cv_box.png***: A box plot showing the cross validation results of meta model.
    - ***vid_YYmmdd_HHMMSS.pkl***: The VID object (for expert usage), timestamped with the current date and time.
    - ***log.txt***: The log message for VID running.

### Parameters ###
__seuratobj_dir__ : str, **requied**
   > The directory of the input rds file (seurat object).

__marker_dir__ : str, **required**
   > The directory of a txt file contains the list of virus biomarkers, with each gene occupying one line.
   > The markers will be applied for the definition of traning set(truely infected and truely uninfected), 
   > while the cells from infected sample with any marker expressed will be considered truely infected.
   > The markers will also be ***excluded*** from modeling.

__clinical_column__ : str, **required**, default = clinical_column
   > The column indicates the sample-level clinical diagnosis of infection, the column should
   > included in the metadata of the input seurat object. 'positive' and 'negative' are two valid values in
   > this column. This column will be applied for training set defination together with 'marker_dir'.
   > Please prepare the column with valid values in metadata and pass the column name when running the code.

__output_dir__ : str, optional, default = './'
   > The output directory, set as current working directory by default.

__feature_dir__ : str, optional, default = None
   > The directory of a text file contains a list of important genes, with each gene occupying one line.
   > If given, ignore feature selection and apply the important genes for modeling, otherwise, the boruta
   > feature selection will be applied to select important genes.

__label_dir__: str, optional, default = None
   > The directory of a text file contains the user-defined label for model training, with each label occupying one line.
   > Three valid labels should be included in the text file: 0 (true negative cell), 1 (true positive cell), and 2 (target cells).
   > If given, ignore the labeling step and apply user-defined label for model construction and prediction. 

__batch_column__ : str, optional, default = None
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

__average__ : str, optional {micro, macro, weighted} , default = weighted 
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

__vidmodel_dir__ : str, optional, default = None
   >  The directory of VID object, perform transfer learning by passing the directory of `vid_YYmmdd_HHMMSS.pkl` from previous output.
   >  The pretained models from the vid object passed will be applied for transfer learning, the output structure is
   >  the same as standard running. Please be cautious that meta model specified should be consistent with the previous running.

### Demo ###

#### NPC-EBV-Lymphocytes ####
create a new work directory for demo:
```bash
mkdir ./demo
```
Download and extract the demo data to the demo directory from dropbox:
```bash
wget --no-check-certificate 'https://www.dropbox.com/scl/fo/3u1ch4939idv6uely16k7/ADsdx0VF-JG2tLQUYYFvhu4?rlkey=eckmr4wpjkiinee6m77fy7rzr&st=hhuhub32&dl=1' -O ./demo.zip
unzip ./demo.zip -d ./demo
```
You can also download from [demo](https://www.dropbox.com/scl/fo/3u1ch4939idv6uely16k7/ADsdx0VF-JG2tLQUYYFvhu4?rlkey=eckmr4wpjkiinee6m77fy7rzr&st=hhuhub32&dl=1) and extract to './demo'.

Running VID:
```bash
run_vid ./demo/data/demo.rds \
--output_dir ./demo \
--marker_dir ./demo/data/EBV_markers.txt \
--clinical_column ebv_status \
--metamodel xgb 
```

#### NPC-EBV-Epithelial ####
Create a new work directory for demo2:
```bash
mkdir ./demo2
```
Download and extract the demo data to the demo directory from dropbox::
```bash
wget --no-check-certificate 'https://www.dropbox.com/scl/fo/1qfrs4izmdxr6pio8a0jx/AMTkwGoQ6samlV7ks3hEu2o?rlkey=gafgpo3j4tg98i0dxun4nf2e1&st=r2s6bcol&dl=1' -O ./demo2.zip
unzip ./demo2.zip -d ./demo2
```
You can also download from [demo2](https://www.dropbox.com/scl/fo/1qfrs4izmdxr6pio8a0jx/AMTkwGoQ6samlV7ks3hEu2o?rlkey=gafgpo3j4tg98i0dxun4nf2e1&st=r2s6bcol&dl=1) and extract to './demo2'.

Running VID:
```bash
run_vid ./demo2/data/demo2.rds \
--output_dir ./demo2 \ 
--marker_dir ./demo2/data/EBV_markers.txt \
--clinical_column EBV_state \
--metamodel mlp
```

### Docker Run ###
Below is the `example` code for docker running:
```
docker run \
-v /path/to/data.rds:/wkdir/input/data.rds \
-v /path/to/output/dir:/wkdir/output \
-v /path/to/marker.txt:/wkdir/input/markers.txt \
-v /path/to/import_genes.txt:/wkdir/input/features.txt \
hwhapply/vid:latest \
--clinical_column clinical_colname \
--optional_argument argument_value \
...
```
Parameter '-v' is applied to map the local directory to container working directory. The usage of '-v' shown below:
```
-v /your/local/dir(file):/container/dir(file)
```
When you execute VID image, you can replace `/your/local/dir(file)` with your local directory, please don't change the `/container/dir(file)`. Modify the container directory will lead to execution failure. 

Specify `clinical_column` after `hwhapply/vid:latest`:
```
docker run
...
hwhapply/vid:latest \
--clinical_column clinical_colname \
...
```
The optional arguments you can specify are listed below:
```
[--batch_column BATCH_COLUMN] [--sample_column SAMPLE_COLUMN] [--test_ratio TEST_RATIO]
[--num_split NUM_SPLIT] [--metamodel METAMODEL] [--threshold THRESHOLD] [--average AVERAGE]
[--random_state RANDOM_STATE] [--n_jobs N_JOBS] [--verbose VERBOSE] [--vidmodel_dir VIDMODEL_DIR]
```

Run VID with docker image on demo data, apply feature (gene) selection and set xgb as meta model:
```bash
docker run \
-v ./demo/data/demo.rds:/wkdir/input/data.rds \
-v ./demo:/wkdir/output \
-v ./demo/data/EBV_markers.txt:/wkdir/input/markers.txt \
hwhapply/vid:latest \
--clinical_column ebv_status \
--metamodel xgb
```
You can also provide the important feature list to skip feature selection:
```bash
docker run \
-v ./demo/data/demo.rds:/wkdir/input/data.rds \
-v ./demo:/wkdir/output \
-v ./demo/data/EBV_markers.txt:/wkdir/input/markers.txt \
-v ./demo/data/important_features.txt:/wkdir/input/features.txt \
hwhapply/vid:latest \
--clinical_column ebv_status \
--metamodel mlp
```
Provide self-defined labels to skip the automatical labeling:
```bash
docker run \
-v ./demo/data/demo.rds:/wkdir/input/data.rds \
-v ./demo:/wkdir/output \
-v ./demo/data/EBV_markers.txt:/wkdir/input/markers.txt \
-v ./demo/data/important_features.txt:/wkdir/input/features.txt \
-v ./demo/data/labels.txt:/wkdir/input/labels.txt \
hwhapply/vid:latest \
--clinical_column ebv_status \
--metamodel mlp
```
The outputs will be saved in the output directory you specified,  in this example the result will be save in `./demo/YYmmdd_HHMMSS` , the structure of docker running output has no different with conda running.

## Transfer VID Model On Unseen Data ##
An object of VID class will be saved as the `vid_YYmmdd_HHMMSS.pkl` in `output` directory of VID output. We can transfer the pre-trained VID model on a new dataset, which can save training time and make the model generalize on unseen data. The pre-trained model can only be applid on the data of same oncovirus and comforms input standard.

### Transfer Learning With Conda ###
To perfrom transer learning in `vid_env` we created, change the input data and specify the directory of pre-trained vid object with argument `vidmodel_dir`.
Perform transfer learning on demo data in `vid_env`:
```bash
run_vid ./demo/data/demo_unseen.rds \
--output_dir ./demo:\
--marker_dir ./demo/data/EBV_markers.txt \
--feature_dir ./demo/data/important_features.txt \
--vidmodel_dir ./demo/data/vid_demo.pkl \
--clinical_column ebv_status \
--metamodel xgb
```

### Transfer Learning With Docker ###
To perform the transfer learning with docker image, exchange the input data with unseen data and map the path of pre-trained vid object to the corresponding directory in the container. Perform transfer learning on the demo dataset with docker image:
```bash
docker run \
-v ./demo/data/demo_unseen.rds:/wkdir/input/data.rds \
-v ./demo:/wkdir/output \
-v ./demo/data/EBV_markers.txt:/wkdir/input/markers.txt \
-v ./demo/data/important_features.txt:/wkdir/input/features.txt \
-v ./demo/data/vid_demo.pkl:/wkdir/input/vid.pkl \
hwhapply/vid:latest \
--clinical_column ebv_status \
--metamodel mlp
```
The output of transfer learning is the same with standard process.

