import pandas as pd
import numpy as np
import time

from sklearn.model_selection import train_test_split 
from sklearn.preprocessing import StandardScaler
import harmonypy as hm
from sklearn.utils.class_weight import compute_class_weight
from sklearn.metrics import make_scorer, confusion_matrix, f1_score, roc_auc_score, accuracy_score, balanced_accuracy_score, precision_score, recall_score
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold, RandomizedSearchCV

import torch 
import torch.nn as nn

from tqdm import tqdm
import anndata
import warnings
warnings.filterwarnings('ignore')  # "error", "ignore", "always", "default", "module" or "once"


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
            if os.path.isfile(self.marker_dir):
                with open(self.marker_dir, 'r') as file:
                    self.markers = [line.strip() for line in file]
            else:
                raise KeyError("Marker file dosen't exist, please double check your marker file directory.")      
        
            # check the if the clinical test column, 'batch' and sample id column included in the metadata table
            essential_columns = [self.clinical_column, self.sample_column]
            for col in essential_columns:
                if col not in self.meta_df.columns:
                    raise KeyError(f'{col} not included in metadata, please add the column.')
                
            # detect the normal control: if a sample only have negative cell, we consider it normal control
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
            print('Labeling...')
            self.markers = [gene for gene in self.markers if gene in self.data_df.columns]
            if len(self.markers) != 0:
                self.meta_df.loc[(self.meta_df[self.clinical_column] == 'positive') &
                                 (self.data_df[self.markers] != 0).any(axis=1), 'label'] = 1
                self.meta_df.loc[self.meta_df[self.sample_column].isin(neg_samples), 'label'] = 0
                self.meta_df['label'].fillna(2, inplace=True)
            else:
                raise KeyError(f'Labeling failed, marker provided not included in the expression data.')
            
            # drop the markers
            self.data_df.drop(self.markers, axis = 1, inplace=True)
            
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
                self.boruta_model.fit(self.X_train.to_numpy(), self.y_train.to_numpy())
                self.features = list(self.X_train.columns[self.boruta_model.support_])
                print(f'Feature selection finished, {len(self.features)} important gene selected.')
            with open(os.path.join(self.output_dir, 'important_genes.txt'), "w") as file:
                for item in self.features:
                    file.write(item + "\n")
            self.X_train = self.X_train[self.features]
            self.X_test = self.X_test[self.features]
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
            model_init_kwargs['xgb']['seed'] = self.random_state
            
            # set class weight for metal model:
            class_weights = compute_class_weight('balanced', classes=np.unique(self.y_train), y=self.y_train)
            class_weights = torch.tensor(class_weights, dtype=torch.float)
            self.pos_weight = class_weights[1] / class_weights[0]
            model_init_kwargs['xgb']['scale_pos_weight'] = float(self.pos_weight)
            model_init_kwargs['mlp']['criterion'] = nn.BCEWithLogitsLoss(pos_weight=self.pos_weight)
                
            
            # set GPU training(CUDA) if available
            self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
            model_init_kwargs['mlp']['device'] = self.device
            
            # set base model and meta model
            if self.pos_weight >=10:
                model_init_kwargs['mlp']['alpha'] = 0.01
                param_grids['param_grid_mlp']['alpha'] = [0.01, 0.05, 0.1, 0.2, 0.5, 1]
            self.models = Model(model_init_kwargs)
            print('Model initialization finished.')
            self.base_models = [self.models.rf, self.models.svc, self.models.knn, self.models.nb, self.models.lgr]
            if self.metamodel == 'mlp':
                self.meta_model = self.models.mlp              
                self.param_grid = param_grids['param_grid_mlp']
            else:
                self.meta_model = self.models.xgb
                self.param_grid = param_grids['param_distributions_xgb']
            
            # define stratified k-fold splitting for hyper-parametering tuning
            self.skf = StratifiedKFold(n_splits = self.num_split, 
                                       shuffle = True,  
                                       random_state = self.random_state)
            
            # define metrics for evaluation
            self.scoring = {
                'accuracy': make_scorer(accuracy_score),
                'balanced_accuracy': make_scorer(balanced_accuracy_score),
                'precision': make_scorer(precision_score, average = self.average, zero_division = 0.0),
                'f1': make_scorer(f1_score, average=self.average),
                'recall': make_scorer(recall_score, average=self.average),
                'roc_auc': make_scorer(roc_auc_score, average=self.average)
                }
                      
        except Exception as e:
            raise ValueError(f'Model initialization failed:\n{e}')
        
    
    def fit(self):
        """
        Fit the stack generalization model on the training set.
        """
        try:
            # *********************************** Step1: Base model training **********************************
            print('Training base models...')
            start_time = time.time()
            # Initialize cv_results dictionary to store the results
            cv_results = {
                'mean_test_accuracy': [],
                'std_test_accuracy': [],
                'mean_test_balanced_accuracy': [],
                'std_test_balanced_accuracy': [],
                'mean_test_recall': [],
                'std_test_recall': [],
                'mean_test_f1': [],
                'std_test_f1': [],
                'mean_test_precision': [],
                'std_test_precision': [],
                'mean_test_roc_auc': [],
                'std_test_roc_auc': []
            }

            # Include split-specific results in cv_results
            for i in range(5):  # Assuming 5 splits
                cv_results[f'split{i}_test_accuracy'] = []
                cv_results[f'split{i}_test_balanced_accuracy'] = []
                cv_results[f'split{i}_test_recall'] = []
                cv_results[f'split{i}_test_f1'] = []
                cv_results[f'split{i}_test_precision'] = []
                cv_results[f'split{i}_test_roc_auc'] = []

            # define variables to save the meta features, scaler 
            meta_features = np.zeros((self.X_train.shape[0], len(self.base_models)))
            test_meta_features = np.zeros((self.X_test.shape[0], len(self.base_models)))
            models = []
            self.scalers = []

            # Track progress of outer loop
            outer_pbar = tqdm(self.skf.split(self.X_train, self.y_train), desc="Training base model", total=self.skf.get_n_splits())

            # Temporary dictionary to hold the metrics for each model and split
            temp_results = {i: {'accuracy': [], 'balanced_accuracy': [], 'recall': [], 'f1': [], 'precision': [], 'roc_auc': []} for i in range(len(self.base_models))}

            for fold_idx, (train_index, val_index) in enumerate(outer_pbar):
                X_train_fold, X_val_fold = self.X_train.iloc[train_index, :], self.X_train.iloc[val_index, :]
                y_train_fold, y_val_fold = self.y_train.iloc[train_index], self.y_train.iloc[val_index]
                columns = X_train_fold.columns
                index_train = X_train_fold.index
                index_val = X_val_fold.index
                
                # perform normalization 
                scaler = StandardScaler()
                X_train_fold = scaler.fit_transform(X_train_fold)
                X_val_fold = scaler.transform(X_val_fold)
                X_test = scaler.transform(self.X_test)
                self.scalers.append(scaler)

                # perform batch correction if batch column is given
                if self.batch_column:
                    ho_train = hm.run_harmony(X_train_fold, self.meta_train.loc[index_train, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                    X_train_fold = pd.DataFrame(ho_train.Z_corr.T, columns=columns, index=index_train)
                    ho_val = hm.run_harmony(X_val_fold, self.meta_train.loc[index_val, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                    X_val_fold = pd.DataFrame(ho_val.Z_corr.T, columns=columns, index=index_val)
                    ho_test = hm.run_harmony(X_test, self.meta_train.loc[self.X_test.index, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                    X_test = pd.DataFrame(ho_test.Z_corr.T, columns=columns, index=self.X_test.index)
                else:
                    X_train_fold = pd.DataFrame(X_train_fold, columns = columns, index = index_train)
                    X_val_fold = pd.DataFrame(X_val_fold, columns = columns, index = index_val)
                    X_test = pd.DataFrame(X_test, columns = columns, index = self.X_test.index)
                    
                models_curr = []

                for i, model in enumerate(self.base_models):
                    # base model training on preprocessed dataset
                    model.fit(X_train_fold, y_train_fold)
                    val_probabilities = model.predict_proba(X_val_fold)[:, 1]  # Get probability of the positive class
                    meta_features[val_index, i] = val_probabilities
                    test_meta_features[:, i] += model.predict_proba(X_test)[:, 1]
                    models_curr.append(model)

                    # Calculate the metrics for the current fold and model
                    val_predictions = model.predict(X_val_fold)
                    accuracy = accuracy_score(y_val_fold, val_predictions)
                    balanced_accuracy = balanced_accuracy_score(y_val_fold, val_predictions)
                    recall = recall_score(y_val_fold, val_predictions, average = self.average)
                    f1 = f1_score(y_val_fold, val_predictions, average = self.average)
                    precision = precision_score(y_val_fold, val_predictions, average = self.average)
                    roc_auc = roc_auc_score(y_val_fold, val_probabilities, average = self.average)
                
                    # Store the results in the temporary dictionary
                    temp_results[i]['accuracy'].append(accuracy)
                    temp_results[i]['balanced_accuracy'].append(balanced_accuracy)
                    temp_results[i]['recall'].append(recall)
                    temp_results[i]['f1'].append(f1)
                    temp_results[i]['precision'].append(precision)
                    temp_results[i]['roc_auc'].append(roc_auc)

                    # Also store the results in the cv_results dictionary for each split
                    cv_results[f'split{fold_idx}_test_accuracy'].append(accuracy)
                    cv_results[f'split{fold_idx}_test_balanced_accuracy'].append(balanced_accuracy)
                    cv_results[f'split{fold_idx}_test_recall'].append(recall)
                    cv_results[f'split{fold_idx}_test_f1'].append(f1)
                    cv_results[f'split{fold_idx}_test_precision'].append(precision)
                    cv_results[f'split{fold_idx}_test_roc_auc'].append(roc_auc)

                models.append(models_curr)

            # Average predictions for the test set 
            self.meta_features = meta_features
            self.test_meta_features = test_meta_features / self.skf.get_n_splits()
            
            # save the trained base models (num_model * num_split)
            self.base_models_trained = models

            # Calculate the mean and std for each metric and model, and store in cv_results
            for i in range(len(self.base_models)):
                cv_results['mean_test_accuracy'].append(np.mean(temp_results[i]['accuracy']))
                cv_results['std_test_accuracy'].append(np.std(temp_results[i]['accuracy']))
                cv_results['mean_test_balanced_accuracy'].append(np.mean(temp_results[i]['balanced_accuracy']))
                cv_results['std_test_balanced_accuracy'].append(np.std(temp_results[i]['balanced_accuracy']))
                cv_results['mean_test_recall'].append(np.mean(temp_results[i]['recall']))
                cv_results['std_test_recall'].append(np.std(temp_results[i]['recall']))
                cv_results['mean_test_f1'].append(np.mean(temp_results[i]['f1']))
                cv_results['std_test_f1'].append(np.std(temp_results[i]['f1']))
                cv_results['mean_test_precision'].append(np.mean(temp_results[i]['precision']))
                cv_results['std_test_precision'].append(np.std(temp_results[i]['precision']))
                cv_results['mean_test_roc_auc'].append(np.mean(temp_results[i]['roc_auc']))
                cv_results['std_test_roc_auc'].append(np.std(temp_results[i]['roc_auc']))

            # Create a DataFrame to summarize the results
            self.basemodel_cv_scores = pd.DataFrame(cv_results, index=['rf', 'svm', 'knn', 'nb', 'lgr'])
            
            print('Base models training finished.')
            
            # *********************************** Step2: Meta Model training **********************************
            print('Training meta model...')
            # prepare the training set for meta_model
            if self.metamodel == 'mlp':
                if not isinstance(self.meta_features, np.ndarray):
                    meta_features = self.meta_features.to_numpy()
                    test_meta_features = self.test_meta_features.to_numpy()
                self.meta_features = torch.tensor(self.meta_features, dtype=torch.float32)
                self.y_train = torch.tensor(self.y_train, dtype=torch.float32).reshape(-1, 1)
                self.test_meta_features = torch.tensor(self.test_meta_features, dtype=torch.float32)
            elif self.metamodel == 'xgb':
                self.meta_features = pd.DataFrame(self.meta_features, columns = ['rf', 'svc', 'knn', 'nb', 'lgr'])
                self.test_meta_features = pd.DataFrame(self.test_meta_features, columns = ['rf', 'svc', 'knn', 'nb', 'lgr'])


            # Train the meta-model using the meta-features
            if self.metamodel == 'mlp':
                grid = GridSearchCV(estimator=self.meta_model,
                                    param_grid=self.param_grid,
                                    cv=self.skf,
                                    scoring=self.scoring,
                                    refit='roc_auc',
                                    verbose=self.verbose,
                                    n_jobs = self.n_jobs
                                    )
            else:
                grid = RandomizedSearchCV(estimator=self.meta_model,
                                          param_distributions=self.param_grid,
                                          cv=self.skf,
                                          scoring=self.scoring,
                                          refit='roc_auc',
                                          verbose=self.verbose,
                                          n_jobs=self.n_jobs,
                                          n_iter=1000  
                                          )

            # hyperparameter tuning 
            self.grid_result = grid.fit(self.meta_features, self.y_train)
            end_time = time.time()
            total_training_time = end_time - start_time
            print(f"Hyperparameter tuning finished, total training time: {total_training_time:.2f} seconds")
            
            # show the feature importance if apply xgbclassifier as metamodel
            if self.metamodel == 'xgb':
                self.xgb_feature_importance()
            
        except Exception as e:
            raise RuntimeError(f'Model training failed:\n{e}')
        
        # perform evaluation 
        try:
            print()
            print('Model Evaluating:')
            self.evaluate()
            print('Model evaluation finished.')
        except Exception as e:
            raise RuntimeError(f'Model Evaluation failed:\n{e}')
        
        # perform prediction
        try:
            print()
            print('Start detection:')
            self.predict()
            print('Detection finished.')
        except Exception as e:
            raise RuntimeError(f'Prediction failed:\n{e}')
        
    def evaluate(self):
        '''
        Evaluate the model on validation set and test set.
        '''
        # *********************************** Base model Evaluation on test set  **********************************
        # Initialize lists to store evaluation scores
        model_names = []
        accuracy_scores = []
        balanced_accuracy_scores = []
        precision_scores = []
        recall_scores = []
        specificity_scores = []
        f1_scores = []
        auc_scores = []
    
        # perform normalization
        scaler = StandardScaler()
        X_train = scaler.fit_transform(self.X_train)
        X_test = scaler.transform(self.X_test)
        
        # perform batch correction
        if self.batch_column:
            ho_train = hm.run_harmony(X_train, self.meta_train.loc[self.X_train.index, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
            X_train = pd.DataFrame(ho_train.Z_corr.T, columns=self.X_train.columns, index=self.X_train.index)
            ho_test = hm.run_harmony(X_test, self.meta_train.loc[self.X_test.index, [self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
            X_test = pd.DataFrame(ho_test.Z_corr.T, columns=self.X_test.columns, index=self.X_test.index)
        else:
            X_train = pd.DataFrame(X_train, columns=self.X_train.columns, index=self.X_train.index)
            X_test = pd.DataFrame(X_test, columns=self.X_test.columns, index=self.X_test.index)

        # Train and evaluate each model
        for model in self.base_models:
            model_name = model.__class__.__name__
            print(f"Training and evaluating {model_name}...")     
            
            # Train the model
            model.fit(X_train, self.y_train)
            
            # Make predictions on the test set
            predictions = model.predict(X_test)
            probabilities = model.predict_proba(X_test)[:, 1]
            
            # Calculate evaluation metrics
            accuracy = accuracy_score(self.y_test, predictions)
            balanced_accuracy = balanced_accuracy_score(self.y_test, predictions)
            precision = precision_score(self.y_test, predictions, average = self.average)
            recall = recall_score(self.y_test, predictions, average = self.average)
            f1 = f1_score(self.y_test, predictions, average = self.average)
            roc_auc = roc_auc_score(self.y_test, probabilities, average = self.average)
            
            # Calculate specificity
            tn, fp, fn, tp = confusion_matrix(self.y_test, predictions).ravel()
            specificity = tn / (tn + fp)
            
            # Append scores to lists
            model_names.append(model_name)
            accuracy_scores.append(accuracy)
            balanced_accuracy_scores.append(balanced_accuracy)
            precision_scores.append(precision)
            recall_scores.append(recall)
            specificity_scores.append(specificity)
            f1_scores.append(f1)
            auc_scores.append(roc_auc)

        # Create a DataFrame to display the scores
        self.basemodel_test_scores = pd.DataFrame({
            'Accuracy': accuracy_scores,
            'Balanced Accuracy': balanced_accuracy_scores,
            'Precision': precision_scores,
            'Recall': recall_scores,
            'Specificity': specificity_scores,
            'F1-Score': f1_scores,
            'AUC': auc_scores
        }, index=model_names)
        
        # *********************************** Meta model evaluation on test set  **********************************
        # Summarize results
        print(self.grid_result.best_estimator_.__class__.__name__)
        print("Best: %f using %s" % (self.grid_result.best_score_, self.grid_result.best_params_))

        # Meta-model predictions on the test set
        final_predictions = self.grid_result.predict(self.test_meta_features)
        final_probabilities = self.grid_result.predict_proba(self.test_meta_features)[:, 1]
        
        # draw and save the confusion matrix
        self.cm_plot(self.y_test, final_predictions)
        
        # draw histogram of the predicted probabilities
        self.histogram(final_probabilities, 'test')
        
        # draw ROC curve of the predicted probabilities
        self.roc_plot(self.y_test, final_probabilities)

        # Calculate evaluation metrics
        accuracy = accuracy_score(self.y_test, final_predictions)
        balanced_accuracy = balanced_accuracy_score(self.y_test, final_predictions)
        precision = precision_score(self.y_test, final_predictions, average = 'weighted')
        recall = recall_score(self.y_test, final_predictions, average = 'weighted')
        f1 = f1_score(self.y_test, final_predictions, average = 'weighted')
        roc_auc = roc_auc_score(self.y_test, final_probabilities, average = 'weighted')

        # Calculate specificity
        tn, fp, fn, tp = confusion_matrix(self.y_test, final_predictions).ravel()
        specificity = tn / (tn + fp)

        # Create a Series for the meta-model's test scores
        meta_test_scores = pd.DataFrame(np.array([accuracy, balanced_accuracy, precision, recall, specificity, f1, roc_auc]).reshape(1,7), 
                                        index=[self.metamodel], 
                                        columns=['Accuracy', 'Balanced Accuracy', 'Precision', 'Recall', 'Specificity', 'F1-Score', 'AUC'])

        # Append the meta-model scores to the scores_df
        self.test_scores = pd.concat([self.basemodel_test_scores, meta_test_scores], axis = 0)
        
        # save the cross validation scores on the test set
        self.test_scores.to_csv(os.path.join(self.output_dir, f'test_scores_{self.average}.csv'))
        
        # *********************************** Cross validation scores on validation set  **********************************
        meta_cv_scores = pd.DataFrame(self.grid_result.cv_results_)
        optimal_meta_cv_score = meta_cv_scores[meta_cv_scores['rank_test_roc_auc'] == 1][0:1][self.basemodel_cv_scores.columns]
        optimal_meta_cv_score.index = [self.metamodel]
        self.val_cv_scores = pd.concat([self.basemodel_cv_scores, optimal_meta_cv_score], axis=0)
        self.val_cv_scores.index = self.test_scores.index
        
        # save the cross validation scores on the test set
        self.val_cv_scores.to_csv(os.path.join(self.output_dir, f'val_cv_scores_{self.average}.csv'))
        
    def predict(self):
        '''
        Predict the infection status of unknown samples.
        '''
        # *********************************** Predict with trained base models  **********************************
        print('Base models prediction...')
        self.unknown_meta_features = np.zeros((self.data_unknown.shape[0], len(self.base_models)))
        i = 0
        for model_fold in tqdm(self.base_models_trained, desc="Predicting:"):
            scaler = self.scalers[i]
            i += 1
            data_unknown = scaler.transform(self.data_unknown[self.features])
            if self.batch_column:
                ho_unknown = hm.run_harmony(data_unknown, self.meta_unknown[[self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
                data_unknown = pd.DataFrame(ho_unknown.Z_corr.T, columns=self.features, index=self.data_unknown.index)
            else:
                data_unknown = pd.DataFrame(data_unknown, columns=self.features, index=self.data_unknown.index)
            for i, model in enumerate(model_fold):
                self.unknown_meta_features[:, i] += model.predict_proba(data_unknown)[:, 1]
        self.unknown_meta_features /= self.skf.get_n_splits()
        self.unknown_meta_features = self.unknown_meta_features.astype('float32')
        print('Base model prediction finished.')

        # *********************************** Predict with trained meta model **********************************
        print('Meta model prediction...')
        pred_proba = list(self.grid_result.predict_proba(self.unknown_meta_features)[:, 1])
        pred = list(self.grid_result.predict(self.unknown_meta_features))
        if self.metamodel == 'mlp':
            pred_proba = [proba[0] for proba in pred_proba]
            pred = [pred[0] for pred in pred]
        print('Meta model prediction finished.')
        
        # concatenate the prediction with meta_data
        df_pred = pd.DataFrame(np.array([pred_proba, pred]).transpose(), index = self.meta_unknown.index, columns = ['infection_probability', 'infection_status'])
        self.meta_unknown = pd.concat([self.meta_unknown, df_pred], axis = 1)
        
        # combine training and unknown metadata
        self.meta_train[['infection_probability']] = self.meta_train[['label']].astype(float)
        self.meta_train[['infection_status']] = self.meta_train[['label']]
        self.meta_train.replace({'infection_status': {1 : 'true_positive', 0 : 'true_negative'}}, inplace=True)
        self.meta_unknown.replace({'infection_status': {1 : 'pred_positive', 0 : 'pred_negative'}}, inplace=True)
        self.meta_df = pd.concat([self.meta_train, self.meta_unknown], axis = 0)
        
        # add another colume with the label generated based on user defined threshold
        if self.threshold:
            self.meta_unknown[f'infection_status_{self.threshold}'] = self.meta_unknown['infection_probability'].apply(lambda x: 'pred_positive' if x >= self.threshold else 'pred_negative')
            self.meta_train[f'infection_status_{self.threshold}'] = self.meta_train[['infection_status']]
        
        # save the updated metadata table
        # if self.h5ad_dir:
        #     self.anndata.obs = self.meta_df
        #     self.anndata.write(self.h5ad_dir)
        # else:
        #     self.meta_df.to_csv(self.meta_dir)
        
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
        self.histogram(self.meta_unknown['infection_probability'], 'unseen')
        
        # visualize the distribution of infection status with umap
        # scaler = StandardScaler()
        # self.data_array_processed = scaler.fit_transform(self.data_df[self.features])
        # if self.batch_column:
        #     ho_all = hm.run_harmony(self.data_array_processed, self.meta_df[[self.batch_column]], vars_use = [self.batch_column], max_iter_harmony=100)
        #     self.data_df_processed = pd.DataFrame(ho_all.Z_corr.T, columns=self.features, index=self.data_df.index)
        # else:
        #     self.data_df_processed = pd.DataFrame(self.data_array_processed, columns=self.features, index=self.data_df.index)
        # self.umap_plot()
    
