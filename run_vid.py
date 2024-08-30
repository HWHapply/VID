import argparse
from VID import *
import pickle
from datetime import datetime

def define_arguments():
	# Create a parser object
	parser = argparse.ArgumentParser(description='Arguments for Virus Infection Detector:')

	# Define optional argument
	parser.add_argument('--h5ad_dir', '-hd', default = None, help = 'The directory of input h5ad file.')
	parser.add_argument('--data_dir', '-dd', default = None, type = str, help = 'The directory of gene expression.')
	parser.add_argument('--meta_dir', '-md', default = None, type = str, help = 'The directory of metadata.')
	parser.add_argument('--output_dir', '-od', default = './output', type = str, help = 'The output directory.')
	parser.add_argument('--marker_dir', '-mkd', default = './markers.txt', type = str, help = 'The markers stores in a txt file(one gene per row).')
	parser.add_argument('--feature_dir', '-fd', default = None, type = str, help = 'The directory of txt file stores the important features(gene).')
	parser.add_argument('--clinical_column', '-cc', default = 'clinical_column', type = str, help = 'The column indicates the infection status in clinical assessment.(Sample level)')
	parser.add_argument('--batch_column', '-bc', default = None, type = str, help = 'The column indicates the batch label that will be used for batch correction(harmony).')
	parser.add_argument('--sample_column', '-sc', default = 'orig.ident', type = str, help = 'The column indicates the sample id.')
	parser.add_argument('--test_ratio', '-tr', default = 0.3, type = float, help = 'The ratio of validating set.')
	parser.add_argument('--num_split', '-ns', default = 5, type = int, help = 'The number of splitting for base model training and hyperparameter tuning for meta model.')
	parser.add_argument('--metamodel', '-mm', default = 'xgb', type = str, help = 'The classifier applied as meta model.')
	parser.add_argument('--threshold', '-threds', default = None, type = float, help = 'The threshold for the decision function of final prediction.')
	parser.add_argument('--average', '-avg', default = 'weighted', type = str, help = 'Define the type of averaging performed on the evaluation scores among different class.')
	parser.add_argument('--random_state', '-rs', default = 42, type = int, help = 'The random state for the reproduction of result.')
	parser.add_argument('--n_jobs', '-threads', default = -1, type = int, help = 'Number of threads applied for parallel excecution.')
	parser.add_argument('--verbose', '-v', default = 2, type = int , help='The verbose mode.')

	# Parse the command-line arguments
	args = parser.parse_args()
	parser.print_help()
	return args


if __name__ == '__main__':
	# Loading arguments
	print("Loading arguments...")
	try:
		args = define_arguments()
		print()
		print("Arguments defined for VID:")
		print()
		args_input = {}
		for arg, value in args.__dict__.items():
			# if value:
			print(f"{arg}: {value}")
			args_input[arg] = value
	except Exception as e:
		raise ValueError("Shutting down due to argument definition error") from e
	vid = VID(args_input)
	vid.fit()
 
	# save the vid model and training time
	current_time = datetime.now().strftime("%Y%m%d_%H%M%S")
	with open(os.path.join(args_input['output_dir'], f'vid_{current_time}.pkl'), "wb") as file:
		pickle.dump(vid, file)
		print('VID model saved!')


