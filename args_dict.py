from Models import MLPClassifier
import torch.optim as optim
import numpy as np
from scipy.stats import uniform, randint
from skorch.callbacks import EarlyStopping


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

#arguments for models initializations
# early stopping initialization for MLP classifier
early_stopping = EarlyStopping(
    monitor='valid_loss',  # Metric to monitor
    patience=5,            # Number of epochs with no improvement after which training will be stopped
    threshold=0.0001,      # Minimum change to qualify as an improvement
    threshold_mode='rel',  # Use relative threshold (i.e., 0.01% improvement)
    lower_is_better=True   # Lower validation loss is better
)

# model initialization arguments
model_init_kwargs = {
    'rf': {
		'n_estimators': 100,
		'criterion': 'gini',
		'max_depth': None,
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
	'svc' : {
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
	'knn' : {
		'n_neighbors' : 5,
		'weights' : 'uniform',
		'algorithm' : 'auto',
		'leaf_size' : 30,
		'p' : 2,
		'metric' : 'minkowski',
		'metric_params' : None,
		'n_jobs' : -1
	},
	'nb': {
		'priors' : None,
		'var_smoothing' : 1e-09
	}, 
     'lgr' : {
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
		'multi_class' : 'auto',
		'verbose' : 0,
		'warm_start' : False,
		'n_jobs' : -1,
		'l1_ratio' : None
	},
 	'xgb' : {# general parameters, 
      	'booster' : 'gbtree',
		'device' : 'cpu',
		'verbosity' : 0,
		'n_jobs' : -1,
		# booster parameters
		'learning_rate': 0.3,
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
  },
  'mlp': {'module':MLPClassifier,
          'optimizer': optim.Adam,
          'max_epochs':100,
          'lr' : 0.01,
          'batch_size': 64,
          'device' : 'cpu',
          'verbose':False,
          'callbacks':[early_stopping],
          'optimizer__weight_decay':0.01
  }
}


# arguments grid for hyperparameter tuning
param_grids = {
    'param_grid_mlp' : {
        'batch_size': [16, 32, 64, 128],
#         'max_epochs': [10, 20, 30, 40, 50, 60, 70, 80, 90, 100],
#         'optimizer': [optim.SGD, optim.RMSprop, optim.Adagrad, optim.Adadelta,
#                       optim.Adam, optim.Adamax, optim.NAdam],
#         'alpha' : [0.01, 0.05, 0.1, 0.2, 0.5, 1],
        'optimizer__lr': np.logspace(-4, -1, 4),
        'optimizer__betas': [(0.85, 0.99)],
#         'optimizer__momentum': np.linspace(0.85, 0.99, 8),
#         'module__weight_init': [init.uniform_, init.normal_, init.zeros_,
#                                 init.xavier_normal_, init.xavier_uniform_,
#                                 init.kaiming_normal_, init.kaiming_uniform_],
#         'module__activation': [nn.Identity, nn.ReLU, nn.ELU, nn.ReLU6,
#                                nn.GELU, nn.Softplus, nn.Softsign, nn.Tanh,
#                                nn.Sigmoid, nn.Hardsigmoid],
#         'module__weight_constraint': [1.0, 2.0, 3.0, 4.0, 5.0],
#         'module__dropout_rate': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
        'module__n_neurons': [16, 32, 48, 64]
    },
    'param_distributions_xgb' : {
        'learning_rate': uniform(0.0001, 0.3),  # Uniform distribution between 0.0001 and 0.3
        'gamma': uniform(0, 5),  # Uniform distribution between 0 and 0.2
        'reg_lambda': uniform(0.1, 10),  # Uniform distribution between 1 and 100
        'reg_alpha': uniform(0, 5),  # Uniform distribution between 0 and 0.1
        'max_depth': randint(3, 6),  # Random integers between 3 (inclusive) and 6 (exclusive)
        'subsample': uniform(0.6, 0.9),  # Uniform distribution between 0.7 and 1.0
        'n_estimators': randint(50, 500)  # Random integers between 10 and 101 (inclusive)
    }
}












