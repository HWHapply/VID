import numpy as np
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from xgboost import XGBClassifier
import torch 
import torch.nn as nn
import torch.nn.init as init
from skorch import NeuralNetClassifier
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
	2. Multi-layer perceptron classifier
"""


# MixupNeuralNetClassifier
class MixupNeuralNetClassifier(NeuralNetClassifier):
    def __init__(self, *args, alpha=1.0, **kwargs):
        super().__init__(*args, **kwargs)
        self.alpha = alpha
        
    def mixup_data(self, x, y, alpha=1.0):
        if alpha > 0:
            lam = np.random.beta(alpha, alpha)
        else:
            lam = 1

        batch_size = x.size()[0]
        index = torch.randperm(batch_size).to(x.device)

        mixed_x = lam * x + (1 - lam) * x[index, :]
        y_a, y_b = y, y[index]
        return mixed_x, y_a, y_b, lam

    def mixup_criterion(self, criterion, pred, y_a, y_b, lam):
        return lam * criterion(pred, y_a) + (1 - lam) * criterion(pred, y_b)

    def get_loss(self, y_pred, y_true, X=None, training=False):
        if not training:
            return super().get_loss(y_pred, y_true, X=X, training=training)
        
        mixed_x, y_a, y_b, lam = self.mixup_data(X, y_true, self.alpha)
        y_pred = self.infer(mixed_x)
        return self.mixup_criterion(self.criterion_, y_pred, y_a, y_b, lam)


# MLP classifier
class MLPClassifier(nn.Module):
    def __init__(self, weight_init=init.kaiming_normal_, 
                 activation=nn.ELU, 
                 dropout_rate=0.5, 
                 weight_constraint=1.0, 
                 n_neurons = 16):
        super().__init__()
        self.layer = nn.Linear(5, n_neurons)
        self.bn = nn.BatchNorm1d(n_neurons)
        self.act = activation()
        self.dropout = nn.Dropout(dropout_rate)
        self.output = nn.Linear(n_neurons, 1)
#        self.prob = nn.Sigmoid()
        # manually init weights
        weight_init(self.layer.weight)
        weight_init(self.output.weight)
        self.weight_constraint = weight_constraint

    def forward(self, x):
        # maxnorm weight before actual forward pass
        with torch.no_grad():
            norm = self.layer.weight.norm(2, dim=0, keepdim=True).clamp(min=self.weight_constraint / 2)
            desired = torch.clamp(norm, max=self.weight_constraint)
            self.layer.weight *= (desired / norm)
        # actual model with skorch
        x = self.bn(self.layer(x))
        x = self.act(x)
        x = self.dropout(x)
        x = self.output(x)
#        x = self.prob(self.output(x))
        return x




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
			if model_name == 'rf':
				return RandomForestClassifier(**kwargs)
			elif model_name == 'svc':
				return SVC(**kwargs)
			elif model_name == 'knn':
				return KNeighborsClassifier(**kwargs)
			elif model_name == 'nb':
				return GaussianNB(**kwargs)
			elif model_name == 'lgr':
				return LogisticRegression(**kwargs)
			elif model_name == 'xgb':
				return XGBClassifier(**kwargs)
			elif model_name == 'mlp':
				if 'alpha' in kwargs:
					return MixupNeuralNetClassifier(**kwargs)
				else:
					return NeuralNetClassifier(**kwargs)
		except Exception as e:
			raise ValueError("Shutting down due to modeling initialization error") from e



	