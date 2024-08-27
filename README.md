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
  * [Simple Usage](#simple-usage)
  * [Parameters](#parameters)
  * [Output](#output)
* [Demo](#demo)
* [Expert Usage](#expert-usage)
  * [16S rRNA gene sequencing data of OSCC patients](#human-microbiome)
  * [Whole metagenomics data of Tara Ocean](#ocean-microbiome)
  * [Transcriptomics data of NASH patients](#human-transcriptome)


## Installation


### Hardware
This model is designed to be executed on standard computational hardware. While the process is feasible on commonly available systems, utilizing more robust computational resources can significantly enhance execution speed and efficiency.

Development and testing were carried out on different OSs:
#### macOS Sonoma Version 14.6.1:
- Processor: 8-core Apple M1 chip
- Memory: 16 GB LPDDR4

#### Ubuntu 22.04:
- Processor: 
- Memory:

#### Windows 10:
- Processor: 
- Memory:

### Software
- R v4.4.1 or newer (https://www.r-project.org)
- Python3 v3.10.14 (https://www.python.org)
- Conda v14.7.1 (https://github.com/conda)

#### R packages
- Seurat (https://satijalab.org/seurat/)
- SeuratDisk (https://github.com/mojaveazure/seurat-disk)

#### python packages
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

#### 3. Configure the environment with a setup script:
```
bash Setup.sh
```

###### If the setup failed, can configure the environment manually with following steps.

Create a conda environment called 'vid_env' for the VID running:
```
conda create -f vid_env.yml -n vid_env
```
Activate the environment:
```
conda activate vid_env
```
Dectivate the VID environment(if required):
```
conda deactivate vid_env
```
Remove the VID environment(if required):
```
conda env remove -n vid_env
```

#### Docker image setup
To provide easier implementation, we provide a Docker image to replace above Equipment setup steps excluding Gephi. Firstly, users should download and install Docker (https://docs.docker.com/engine/install/) and then setup the xMarkerFinder computational environment. All scripts in the Procedure part below should be executed within the Docker container created from the xMarkerFinder Docker image.

```
$ docker pull tjcadd2022/xmarkerfinder:1.0.16
$ docker run -it -v $(pwd):/work tjcadd2022/xmarkerfinder:1.0.16 /bin/bash  
```
```
-it Run containers in an interactive mode, allowing users to execute commands and access files within the docker container.  
-v Mounts a volume between present working directory in your local machine to the /work directory in the docker container.  
```

## User tutorial
### Simple Usage ###
```

``` 
- Input files:
Seurat object: merged microbial count profile of all datasets.  
- Output files:  
normalized_abundance.txt: normalized abundance profile of input dataset. Normalized abundance profiles are used as input files for all subsequent analyses, except for Step 11, which requires raw count file.

### Parameters ###

__estimator__ : object
   > A supervised learning estimator, with a 'fit' method that returns the
   > feature_importances_ attribute. Important features must correspond to
   > high absolute values in the feature_importances_.

__n_estimators__ : int or string, default = 1000
   > If int sets the number of estimators in the chosen ensemble method.
   > If 'auto' this is determined automatically based on the size of the
   > dataset. The other parameters of the used estimators need to be set
   > with initialisation.

__perc__ : int, default = 100
   > Instead of the max we use the percentile defined by the user, to pick
   > our threshold for comparison between shadow and real features. The max
   > tends to be too stringent. This provides a finer control over this. The
   > lower perc is the more false positives will be picked as relevant but
   > also the less relevant features will be left out. The usual trade-off.
   > The default is essentially the vanilla Boruta corresponding to the max.

__alpha__ : float, default = 0.05
   > Level at which the corrected p-values will get rejected in both correction
   steps.

__two_step__ : Boolean, default = True
  > If you want to use the original implementation of Boruta with Bonferroni
  > correction only set this to False.

__max_iter__ : int, default = 100
   > The number of maximum iterations to perform.

__verbose__ : int, default=0
   > Controls verbosity of output.

### Output ###
20240827_174623
├── data
│   └── data.h5ad
└── output
    ├── Confusion_Matrix_mlp.png
    ├── ROC_Curve_mlp.png
    ├── important_genes.txt
    ├── pred_proba_hist_mlp_test.png
    ├── pred_proba_hist_mlp_unseen.png
    ├── test_scores_weighted.csv
    ├── val_cv_scores_weighted.csv
    └── vid_240827_1748.pkl

## Demo ##


### Seurat Object Version < V5 ###


### Seurat Object V5 ###


### MLP as Meta model ###

## Expert Usage ##
