import numpy as np
from scipy.stats import uniform, randint



# feature selection arguments
args_fs = {
    'estimator' : {
        'n_estimators' : 100, 
        'criterion' : 'gini', 
		'max_depth' : None, 
		'min_samples_split' : 2, 
		'min_samples_leaf' : 1, 
		'min_weight_fraction_leaf' : 0.0, 
		'max_features' : 'sqrt', 
		'max_leaf_nodes' : None, 
		'min_impurity_decrease' : 0.0, 
		'bootstrap' : True, 
		'oob_score' : False, 
		'n_jobs' : -1, 
		'random_state' : 42, 
		'verbose' : 0, 
		'warm_start' : False, 
		'class_weight' : 'balanced', 
		'ccp_alpha' : 0.0, 
		'max_samples' : None
    },
	'boruta' : {
		'estimator' : None,
		'perc' : 100,
		'alpha' : 0.03,
		'two_step' : True,
		'max_iter' : 100,
		'verbose' : 0,
		'random_state' : 42
	}
}

# model initialization arguments
model_init_kwargs = {
    'RF': {
		'n_estimators': 100,
		'criterion': 'gini',
		'max_depth': 5,
		'min_samples_split': 2,
		'min_samples_leaf': 1,
		'min_weight_fraction_leaf': 0.0,
		'max_features': 'sqrt',
		'max_leaf_nodes': None,
		'min_impurity_decrease': 0.0,
		'bootstrap': True,
		'oob_score': False,
		'n_jobs': -1,
		'random_state': 42,
		'verbose': 0,
		'warm_start': False,
		'class_weight': 'balanced',
		'ccp_alpha': 0.0,
		'max_samples': None
	},
	'SVM' : {
		'C' : 1.0,
		'kernel' : 'rbf',
		'degree' : 3,
		'gamma' : 'scale',
		'coef0' : 0.0,
		'shrinking' : True,
		'probability' : True,
		'tol' : 0.001,
		'cache_size' : 200,
		'class_weight' : 'balanced',
		'verbose' : False,
		'max_iter' : -1,
		'decision_function_shape' : 'ovr',
		'break_ties' : False,
		'random_state' : 42
	},
	'KNN' : {
		'n_neighbors' : 5,
		'weights' : 'uniform',
		'algorithm' : 'auto',
		'leaf_size' : 30,
		'p' : 2,
		'metric' : 'minkowski',
		'metric_params' : None,
		'n_jobs' : -1
	},
	'GNB': {
		'priors' : None,
		'var_smoothing' : 1e-09
	}, 
     'LGR' : {
		'penalty' : 'l2',
		'dual' : False,
		'tol' : 0.0001,
		'C' : 1.0,
		'fit_intercept': True,
		'intercept_scaling' : 1,
		'class_weight' : 'balanced',
		'random_state' : 42,
		'solver' : 'lbfgs',
		'max_iter' : 5000,
		'verbose' : 0,
		'warm_start' : False,
		'n_jobs' : -1,
		'l1_ratio' : None
	},
 	'XGB' : {# general parameters, 
      	'booster' : 'gbtree',
		'device' : 'cpu',
		'verbosity' : 0,
		'n_jobs' : -1,
		# booster parameters
		'learning_rate': 0.1,
		'gamma' : 0,
		'max_depth': 5,
		'alpha': 0,
		'min_child_weight' : 1,
		'max_delta_step' : 0,
		'subsample' : 1,
		'sampling_method' : 'uniform',
		'colsample_bytree' : 1, 
		'colsample_bylevel' : 1, 
		'colsample_bynode' : 1,
		'lambda' : 1,
		'tree_method' : 'auto',
		'scale_pos_weight' : 1,
		'n_estimators' : 100,
		# learning task parameter
		'objective':'binary:logistic',
		'eval_metric' : 'auc',
		'seed' : 42
  }
}


# arguments grid for hyperparameter tuning
param_grids = {
    # Base: RandomForest
    'RF__n_estimators': randint(50, 200),
    'RF__criterion': ['gini', 'entropy', 'log_loss'],
    'RF__max_depth': randint(3, 10),
    'RF__min_samples_split': randint(2, 10),
    'RF__min_samples_leaf': randint(1, 4),

    # Base: SVC (in pipeline with StandardScaler)
    'SVM__C': uniform(0.1, 10),
    'SVM__gamma': ['scale', 'auto'],
    'SVM__kernel': ['poly', 'rbf', 'sigmoid'],

    # Base: KNN
    'KNN__n_neighbors': randint(3, 10),
    'KNN__weights': ['uniform', 'distance'],
    'KNN__p': [1, 2],

    # Base: GaussianNB — only one tunable param
    'GNB__var_smoothing': uniform(1e-11, 1e-8),

    # Base: Logistic Regression
    'LGR__C': uniform(0.01, 10),
    
    # Final: XGB
	'final_estimator__learning_rate': uniform(0.0001, 0.3),  # Uniform distribution between 0.0001 and 0.3
	'final_estimator__gamma': uniform(0, 5),  # Uniform distribution between 0 and 0.2
	'final_estimator__reg_lambda': uniform(0.1, 10),  # Uniform distribution between 1 and 100
	'final_estimator__reg_alpha': uniform(0, 5),  # Uniform distribution between 0 and 0.1
	'final_estimator__max_depth': randint(3, 6),  # Random integers between 3 (inclusive) and 6 (exclusive)
	'final_estimator__subsample': uniform(0.6, 0.3),  # Uniform distribution between 0.7 and 1.0
	'final_estimator__n_estimators': randint(50, 200)  # Random integers between 10 and 101 (inclusive)
}












