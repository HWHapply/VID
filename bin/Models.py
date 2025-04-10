import numpy as np
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
import warnings
warnings.filterwarnings("ignore")


"""
Models:
Base models:
	1. Random forest
	2. Support vector machine
	3. K nearest neighbor
	4. Naive bayes
	5. Logistic regression
 Meta models:
	1. Extreme gradient boosting
"""

class Model:
	"""
    This class initialize the ML models with given arguments.

    Arguments:
        args_grid(dict{str(model): dict{str(arg1): .., str(arg2): ..}, ..}) : The initial arguments.
    """

	def __init__(self, kwargs_grid):
		for key, kwargs in kwargs_grid.items():
			model = self.init_model(key, kwargs)
			print(f'Model {key} initialized.')
			setattr(self, key, model)

	def init_model(self, model_name, kwargs):
		"""
        Initialize the model with the give arguments.

        Arguments:
            model_name(str): the model name.
            kwargs(dict): the initial arguments.

        Return:
            model(sklearn.classifier): Initialized model.
        """
		try:
			if model_name == 'RF':
				return RandomForestClassifier(**kwargs)
			elif model_name == 'SVM':
				return SVC(**kwargs)
			elif model_name == 'KNN':
				return KNeighborsClassifier(**kwargs)
			elif model_name == 'GNB':
				return GaussianNB(**kwargs)
			elif model_name == 'LGR':
				return LogisticRegression(**kwargs)
			elif model_name == 'XGB':
				return XGBClassifier(**kwargs)

		except Exception as e:
			raise ValueError("Shutting down due to modeling initialization error") from e



	