import argparse
from VID import *
import pickle
from datetime import datetime


def define_arguments():
	# Create a parser object
	parser = argparse.ArgumentParser()

	# Define optional argument
	parser.add_argument('--h5ad_dir', '-hd', default = None, help = 'The directory of input h5ad file.')
	parser.add_argument('--data_dir', '-dd', default = None, type = str, help = 'The directory of gene expression.')
	parser.add_argument('--meta_dir', '-md', default = None, type = str, help = 'The directory of metadata.')
	parser.add_argument('--output_dir', '-od', default = './output', type = str, help = 'The output directory.')
	parser.add_argument('--n_iter', '-nit', default = 100, type = int, help = 'Number of iteration applied in randomsearchcv.')
	parser.add_argument('--marker_dir', '-mkd', default = None, type = str, help = 'The markers stores in a txt file(one gene per row).')
	parser.add_argument('--feature_dir', '-fd', default = None, type = str, help = 'The directory of txt file stores the important features(gene).')
	parser.add_argument('--label_dir', '-ld', default = None, type = str, help = 'The directory of txt file stores the pre-defined labels.')
	parser.add_argument('--clinical_column', '-cc', default = 'clinical_column', type = str, help = 'The column indicates the infection status in clinical assessment.(Sample level)')
	parser.add_argument('--batch_column', '-bc', default = None, type = str, help = 'The column indicates the batch label that will be used for batch correction(harmony).')
	parser.add_argument('--sample_column', '-sc', default = 'orig.ident', type = str, help = 'The column indicates the sample id.')
	parser.add_argument('--test_ratio', '-tr', default = 0.3, type = float, help = 'The ratio of validating set.')
	parser.add_argument('--num_split', '-ns', default = 5, type = int, help = 'The number of splitting for base model training and hyperparameter tuning for meta model.')
	parser.add_argument('--threshold', '-threds', default = None, type = float, help = 'The threshold for the decision function of final prediction.')
	parser.add_argument('--random_state', '-rs', default = 42, type = int, help = 'The random state for the reproduction of result.')
	parser.add_argument('--n_jobs', '-threads', default = -1, type = int, help = 'Number of threads applied for parallel excecution.')
	parser.add_argument('--verbose', '-v', default = 1, type = int , help='The verbose mode.')

	# Parse the command-line arguments
	args = parser.parse_args()
	#parser.print_help()
	return args


if __name__ == '__main__':
	# Loading arguments
	print("Loading arguments...")
	try:
		args = define_arguments()
		print()
		print("Arguments defined for VID:")
		args_input = {}
		for arg, value in args.__dict__.items():
			# if value:
			print(f"{arg}: {value}")
			args_input[arg] = value
	except Exception as e:
		raise ValueError("Shutting down due to argument definition error") from e

	# initialize and train VID
	vid = VID(args_input)
	vid.fit()
  
	# perform evaluation 
	try:
		print()
		print('Model Evaluating:')
		vid.evaluate()
		print('Model evaluation finished.')
	except Exception as e:
		raise RuntimeError(f'Model Evaluation failed:\n{e}')
	
	# perform prediction
	try:
		print()
		print('Start detection...')
		vid.predict_unknown()
		print('Detection finished.')
	except Exception as e:
		raise RuntimeError(f'Prediction failed:\n{e}')
 


