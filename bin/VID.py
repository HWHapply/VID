import pandas as pd
import numpy as np
import time

from sklearn.model_selection import train_test_split 
from sklearn.preprocessing import StandardScaler
import harmonypy as hm
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import make_scorer, confusion_matrix, f1_score, roc_auc_score, accuracy_score, balanced_accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, RandomizedSearchCV
from sklearn.utils.validation import check_is_fitted
from sklearn.ensemble import StackingClassifier

import anndata
import warnings
warnings.filterwarnings('ignore')  



from Utils import Utils_Model
from sklearn.preprocessing import StandardScaler
from Models import Model
import os
from args_dict import args_fs, model_init_kwargs,  param_grids

class VID(Utils_Model):
    def __init__(self, kwargs):
        """
        Construct Simple machine learning pipeline.

        Args:
            kwargs (dict): Dictionary contains training configuration.
        """
        # define instance attribute
        for key, value in kwargs.items():
            setattr(self, key, value)
            
        # prepare output directory
        if not os.path.isdir(self.output_dir):
            os.makedirs(self.output_dir)
                      
        # dataset preparation
        self.preprocessing()

        # initialize model with specified arguments
        self.model_init()
    
        
    def preprocessing(self):
        try:
            # load the dataset
            print('Loading dataset...')
            # load h5ad file if anndata path provided
            if self.h5ad_dir:
                self.anndata = anndata.read_h5ad(self.h5ad_dir)
                self.data_df = pd.DataFrame(data=self.anndata.raw[:, self.anndata.var_names].X.toarray(), index = self.anndata.raw.obs_names, columns = self.anndata.var_names)
                self.meta_df = self.anndata.obs
            # load gene expression and metadata if anndata path not provided
            elif os.path.isfile(self.data_dir) and os.path.isfile(self.meta_dir):
                self.data_df = pd.read_csv(self.data_dir, index_col=0)
                self.meta_df = pd.read_csv(self.meta_dir, index_col=0)
                if self.data_df.shape[0] != self.meta_df.shape[0]:
                    if self.data_df.transpose().shape[0] == self.meta_df.shape[0]:
                        self.data_df = self.data_df.transpose()
                    else:
                        raise ImportError("The gene expression dosen't match meta data.")
            print('Dataset loading finished.')
            
            # load the marker list of the virus
            if self.marker_dir:
                with open(self.marker_dir, 'r') as file:
                    self.markers = [line.strip() for line in file]
                    self.markers = [gene for gene in self.markers if gene in self.data_df.columns]
        
            # check the if the clinical test column and sample id column included in the metadata table
            essential_columns = [self.clinical_column, self.sample_column]
            for col in essential_columns:
                if col not in self.meta_df.columns:
                    raise KeyError(f'{col} not included in metadata, please add the column.')
                
            # convert the clinical test value to lowercase and check if only 'positive' and 'negative' included
            self.meta_df[self.clinical_column] = self.meta_df[self.clinical_column].str.lower()
            if not self.meta_df[self.clinical_column].isin(['positive', 'negative']).all():
                raise ValueError("Only 'positive' and 'negative' allowed in the clinical column, unexpected value provided.")
    
                
            # detect the normal control: if a sample only have negative cell, we consider it normal control
            if not self.label_dir:
                neg_value_counts = self.meta_df.loc[self.meta_df[self.clinical_column] == 'negative', self.sample_column].value_counts()
                neg_samples = []
                for sid in neg_value_counts.index:
                    if neg_value_counts.loc[sid] == len(self.meta_df[self.meta_df[self.sample_column] == sid]):
                        neg_samples.append(sid)
                if len(neg_samples) == 0:
                    raise ValueError(f'Cannot detect the normal control(true negative sample) in provided dataset.')
                else:
                    print(f'{len(neg_samples)} normal controls detected.')
            
            # labeling based on the clinical test and gene expression
            # true positive : clinical positive + at least one marker expressed
            # true negative : clinical negative 
            # unknown sample : clinical positive + no marker expressed
            if self.label_dir:
                print('Applying the pre-defined labels.')
                with open(self.label_dir, 'r') as file:
                    self.labels = [int(line.strip()) for line in file]
                self.meta_df['label'] = self.labels
                file.close()
            else:
                print('Labeling...')
                if len(self.markers) != 0:
                    self.meta_df.loc[(self.meta_df[self.clinical_column] == 'positive') &
                                    (self.data_df[self.markers] != 0).any(axis=1), 'label'] = 1
                    self.meta_df.loc[self.meta_df[self.sample_column].isin(neg_samples), 'label'] = 0
                    self.meta_df['label'].fillna(2, inplace=True)
                else:
                    raise KeyError(f'Labeling failed, marker provided not included in the expression data.')
            
            # drop the markers
            # if self.marker_dir:
            #     self.data_df.drop(self.markers, axis = 1, inplace=True)
            
            self.data_train = self.data_df[self.meta_df['label'] != 2]
            self.meta_train = self.meta_df.loc[self.data_train.index, :]
            self.data_unknown = self.data_df[self.meta_df['label'] == 2]
            self.meta_unknown = self.meta_df.loc[self.data_unknown.index, :]
            print('Dataset labeling finished.')
            
            # perform stratified spliting on training set
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
                self.data_train,
                self.meta_train['label'],
                test_size=self.test_ratio,
                random_state=self.random_state,
                shuffle=True,
                stratify=self.meta_train['label'])
            print(f'Dataset stratified splitting finished:')
            print(f'Class distribution among {self.X_train.shape[0]} training samples:')
            print(self.y_train.value_counts())
            print(f'Class distribution among {self.X_test.shape[0]} testing samples:')
            print(self.y_test.value_counts())

            # perform normalization
            scaler = StandardScaler()
            self.X_train_norm = scaler.fit_transform(self.X_train)
            self.X_test_norm = scaler.transform(self.X_test)
            self.data_unknown_norm = scaler.transform(self.data_unknown)

            # perform batch correction if specified
            if self.batch_column:
                ho_train = hm.run_harmony(self.X_train_norm, self.meta_train.loc[self.X_train.index, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                self.X_train_norm = pd.DataFrame(ho_train.Z_corr.T, columns=self.X_train.columns, index=self.X_train.index)
                ho_test = hm.run_harmony(self.X_test_norm, self.meta_train.loc[self.X_test.index, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                self.X_test_norm = pd.DataFrame(ho_test.Z_corr.T, columns=self.X_test.columns, index=self.X_test.index)
                ho_unknown = hm.run_harmony(self.data_unknown_norm, self.meta_unknown[[self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                self.data_unknown_norm = pd.DataFrame(ho_unknown.Z_corr.T, columns=self.data_unknown.columns, index=self.data_unknown.index)
            else:
                self.X_train_norm = pd.DataFrame(self.X_train_norm, columns=self.X_train.columns, index=self.X_train.index)
                self.X_test_norm = pd.DataFrame(self.X_test_norm, columns=self.X_test.columns, index=self.X_test.index)
                self.data_unknown_norm = pd.DataFrame(self.data_unknown_norm, columns=self.data_unknown.columns, index=self.data_unknown.index)
            
            # feature selection
            if self.feature_dir:
                print('Loading important features...')
                with open(self.feature_dir, 'r') as file:
                    self.features = [line.strip() for line in file]
                print('Important features loaded.')
            else:
                print('Feature selecting...')
                args_fs['estimator']['random_state'] = self.random_state
                args_fs['boruta']['random_state'] = self.random_state
                self.boruta_model = self.boruta()
                self.boruta_model.fit(self.X_train_norm.to_numpy(), self.y_train.to_numpy())
                self.features = list(self.X_train.columns[self.boruta_model.support_])
                print(f'Feature selection finished, {len(self.features)} important gene selected.')
            with open(os.path.join(self.output_dir, 'important_genes.txt'), "w") as file:
                for item in self.features:
                    file.write(item + "\n")
            
                                
        except Exception as e:
            raise RuntimeError(f'Dataset preparation failed: {e}')


    def model_init(self):
        try: 
            # set random states for all stochastic model
            # Loop through features selection arguments and update 'random_state' and 'n_jobs'
            for key in args_fs:
                if 'random_state' in args_fs[key]:
                    args_fs[key]['random_state'] = self.random_state
                    
            for key in model_init_kwargs:
                if 'random_state' in model_init_kwargs[key]:
                    model_init_kwargs[key]['random_state'] = self.random_state
                if 'n_jobs' in model_init_kwargs[key]:
                    model_init_kwargs[key]['n_jobs'] = self.n_jobs
            
            # set random state for xgb
            model_init_kwargs['XGB']['seed'] = self.random_state
            model_init_kwargs['XGB']['verbosity'] = self.verbose
            
            # set class weight for metal model:
            class_weights = compute_class_weight('balanced', classes=np.unique(self.y_train), y=self.y_train)
            self.pos_weight = class_weights[1] / class_weights[0]
            model_init_kwargs['XGB']['scale_pos_weight'] = float(self.pos_weight)
            
            # initialize base model and meta model
            self.models = Model(model_init_kwargs)
            print('Model initialization finished.')
            self.base_models = [self.models.RF, self.models.SVM, self.models.KNN, self.models.GNB, self.models.LGR]
            self.meta_model = self.models.XGB
            
            # define stratified k-fold splitting for hyper-parametering tuning
            self.skf = StratifiedKFold(n_splits = self.num_split, 
                                       shuffle = True,  
                                       random_state = self.random_state)
                      
        except Exception as e:
            raise ValueError(f'Model initialization failed:\n{e}')
        
    
    def fit(self):
        """
        Fit the stack generalization model on the training set.
        """
        try:
            print('')
            print('Start model training...')
            start_time = time.time()
            
            # initialize stack generalization model
            self.base_estimators = [
                ('RF', self.base_models[0]),
                ('SVM', self.base_models[1]),
                ('KNN', self.base_models[2]),
                ('GNB', self.base_models[3]),
                ('LGR', self.base_models[4])
            ]
            
            stack = StackingClassifier(
                estimators=self.base_estimators,
                final_estimator=self.meta_model,
                cv=self.skf,
                stack_method='auto',
                passthrough=False, 
                n_jobs=self.n_jobs,
                verbose=self.verbose,
            )

            # Train the model and apply hyperparameter tunning
            grid = RandomizedSearchCV(estimator=stack,
                                      param_distributions=param_grids,
                                      cv=self.skf,
                                      refit='roc_auc',
                                      verbose=self.verbose,
                                      n_jobs=self.n_jobs,
                                      n_iter=self.n_iter,
                                      random_state=self.random_state
                                      )

            self.grid_result = grid.fit(self.X_train_norm[self.features], self.y_train)
            end_time = time.time()
            total_training_time = end_time - start_time
            print(f"Model training finished, total training time: {total_training_time:.2f} seconds")
            
            # show the feature importance of xgb 
            self.xgb_feature_importance()
            
        except Exception as e:
            raise RuntimeError(f'Model training failed:\n{e}')
        
      
    def evaluate(self):
        '''
        Evaluate the model on validation set and test set.
        '''
        # check the model fitting status
        
        # *********************************** Base model Evaluation on test set  **********************************
        # Initialize lists to store evaluation scores
        self.pred_proba = {}
        self.clf_list = []
    
        print('Start model evaluation...')
        # Train and evaluate each model
        for name, model in self.base_estimators:
            model.fit(self.X_train_norm, self.y_train)
            y_prob = model.predict_proba(self.X_test_norm)[:, 1]      
            # save the predicted probabilities of base models
            self.pred_proba[name] = y_prob
            self.clf_list.append((model, name))
        
        # *********************************** VID evaluation on test set  **********************************
        stack_pred_proba = self.grid_result.best_estimator_.predict_proba(self.X_test_norm[self.features])[:, 1]
        stack_pred = self.grid_result.best_estimator_.predict(self.X_test_norm[self.features])
        
        # save the predicted probability of meta model
        self.pred_proba['VID'] = stack_pred_proba
        self.clf_list.append((self.grid_result.best_estimator_, 'VID'))
        
        print('Drawing evaluation plots...')
        # draw and save the confusion matrix
        self.cm_plot(self.y_test, stack_pred)
        
        # draw ROC and precision-recall curve of the predicted probabilities
        self.roc_pr_plot(self.y_test, stack_pred_proba)
        
        # draw the forest plot
        self.ci_dict = self.evaluate_models_with_bootstrap()
        self.forest_plot()
        self.combine_images_pillow()
        
        # draw the calibration plot
        #self.calibration_plot()
        
        print('Model Evaluation finished.')
        
        
    def predict_unknown(self):
        '''
        Predict the infection status of unknown samples.
        '''        
        # *********************************** Predict optimal Model  **********************************

        stack_pred_proba = self.grid_result.best_estimator_.predict_proba(self.data_unknown_norm[self.features])[:, 1]
        stack_pred = self.grid_result.best_estimator_.predict(self.data_unknown_norm[self.features])
        
        # concatenate the prediction with meta_data
        df_pred = pd.DataFrame(np.array([stack_pred_proba, stack_pred]).transpose(), index = self.meta_unknown.index, columns = ['infection_probability', 'infection_status'])
        self.meta_unknown = pd.concat([self.meta_unknown, df_pred], axis = 1)
        
        # combine training and unknown metadata
        self.meta_train[['infection_probability']] = self.meta_train[['label']].astype(float)
        self.meta_train[['infection_status']] = self.meta_train[['label']]
        self.meta_train.replace({'infection_status': {1 : 'true_positive', 0 : 'true_negative'}}, inplace=True)
        self.meta_unknown.replace({'infection_status': {1 : 'pred_positive', 0 : 'pred_negative'}}, inplace=True)
        # add another colume with the label generated based on user defined threshold
        if self.threshold:
            self.meta_unknown[f'infection_status_{self.threshold}'] = self.meta_unknown['infection_probability'].apply(lambda x: 'pred_positive' if x >= self.threshold else 'pred_negative')
            self.meta_train[f'infection_status_{self.threshold}'] = self.meta_train[['infection_status']]
        self.meta_df = pd.concat([self.meta_train, self.meta_unknown], axis = 0)
        
        
        # save the metadata table
        output_data_dir = os.path.join(os.path.dirname(self.output_dir), 'data')
        if os.path.isfile(os.path.join(output_data_dir, 'dmatrix.csv')):
            os.remove(os.path.join(output_data_dir, 'dmatrix.csv'))
        if os.path.isfile(os.path.join(output_data_dir, 'data.h5ad')):
            os.remove(os.path.join(output_data_dir, 'data.h5ad'))
        if not self.meta_dir:
            self.meta_dir = os.path.join(output_data_dir, 'metadata.csv')
        self.meta_df.to_csv(self.meta_dir)
        print('Metatable saved!')

        # draw histogram of the predicted probabilities
        self.histogram()

        
 

     


