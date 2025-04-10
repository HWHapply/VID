import matplotlib.pyplot as plt
import matplotlib
from sklearn.metrics import confusion_matrix, roc_auc_score, ConfusionMatrixDisplay, roc_curve, RocCurveDisplay, PrecisionRecallDisplay, precision_score, recall_score, f1_score, balanced_accuracy_score
from boruta import BorutaPy
from args_dict import args_fs
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import resample

import os
import xgboost as xgb
import seaborn as sns
from zepid.graphics import EffectMeasurePlot
from matplotlib.gridspec import GridSpec
from sklearn.calibration import CalibratedClassifierCV, CalibrationDisplay
from PIL import Image
import numpy as np
import pandas as pd
np.int = np.int32
np.float = np.float64
np.bool = np.bool_


"""
Utils_model:
1. Feature selection
2. Confusion Matrix plot
3. ROC Curve plot
4. Histogram of predicted probabilities
5. XGB feature importance bar plot
6. Umap projection
"""


class Utils_Model:
    '''
    Prepare util functions.
    '''
    def boruta(self):
        '''
        Initialize and return a boruta feature selector 
        '''
        self.fs_estimator = RandomForestClassifier(**args_fs['estimator'])
        args_fs['boruta']['estimator'] = self.fs_estimator
        args_fs['boruta']['verbose'] = self.verbose
        boruta_selector = BorutaPy(**args_fs['boruta'])
        return boruta_selector
    
    def cm_plot(self, label, pred):
        """
        Draw the confusion matrix.
        
        Parameters:
        - pred: predictions of model on testing set
        - label: true class label of testing set

        Returns:
        - A confusion matrix plot with larger values.
        """
        cm = confusion_matrix(label, pred)

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(5, 4))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['Uninfected', 'Infected'])

        # Plot confusion matrix
        disp.plot(include_values=True, cmap='Blues', values_format='d', ax=ax)


        # Remove plot spines (grid box lines)
        for spine in ax.spines.values():
            spine.set_visible(False)
            

        # Set tick label fonts and sizes
        ax.set_xticklabels(disp.display_labels, fontsize=8, fontname='Arial', va='center')
        ax.set_yticklabels(disp.display_labels, fontsize=8, fontname='Arial', rotation = 90, va='center')

        # Keep the tick labels 
        ax.tick_params(axis='x', which='both', length=0, pad = 8)  # Removes x-axis ticks
        ax.tick_params(axis='y', which='both', length=0)  # Removes y-axis ticks
        
        # Set axis labels
        ax.set_xlabel('Prediction', fontsize=8, fontname='Arial', loc='center')
        ax.set_ylabel('Ground Truth', fontsize=8, fontname='Arial', loc='center')

        # Customize number annotations inside matrix
        for text in disp.text_.ravel():
            text.set_fontsize(8)
            text.set_fontname('Arial')

        # Colorbar font customization
        cbar = disp.im_.colorbar
        cbar.ax.tick_params(labelsize=8)
        for label in cbar.ax.get_yticklabels():
            label.set_fontname('Arial')
            
        # Remove the colorbar border
        for spine in cbar.ax.spines.values():
            spine.set_visible(False)

        plt.tight_layout()
        # Save figure
        output_path = os.path.join(self.output_dir, 'Confusion_matrix_test.png')
        plt.savefig(output_path, bbox_inches="tight", dpi=600)
        plt.close()

    def roc_pr_plot(self, y_true, y_pred_proba):
        """
        Draw combined ROC and Precision-Recall (P-R) plots side by side.

        Args:
            y_true : np.array, the list of true labels.
            y_pred_proba : np.array, the list of predicted probabilities.
        """
        # Set global font
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 8

        # Create figure with 2 subplots
        fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=600)

        # === ROC Curve ===
        roc_display = RocCurveDisplay.from_predictions(
            y_true,
            y_pred_proba,
            name="VID",
            color="blue",
            ax=axes[0],
            plot_chance_level=True
        )

        axes[0].set_xlabel("False Positive Rate", fontsize=8, fontname='Arial')
        axes[0].set_ylabel("True Positive Rate", fontsize=8, fontname='Arial')
        axes[0].set_title("ROC Curve", fontsize = 8, fontname = 'Arial')
        axes[0].set_xlim(-0.05, 1.05)
        axes[0].set_ylim(-0.05, 1.05)

        for label in axes[0].get_xticklabels() + axes[0].get_yticklabels():
            label.set_fontname('Arial')
            label.set_fontsize(8)

        roc_legend = axes[0].get_legend()
        if roc_legend:
            roc_legend.set_title(None)
            for text in roc_legend.get_texts():
                text.set_fontname('Arial')
                text.set_fontsize(8)

        # === PR Curve ===
        pr_display = PrecisionRecallDisplay.from_predictions(
            y_true,
            y_pred_proba,
            name="VID",
            color="blue",
            ax=axes[1]
        )

        baseline = np.sum(y_true) / len(y_true)
        axes[1].hlines(baseline, 0, 1, color='gray', linestyle='--', linewidth=1,
                    label=f'Chance level (AP = {baseline:.2f})')

        axes[1].set_xlabel("Recall", fontsize=8, fontname='Arial')
        axes[1].set_ylabel("Precision", fontsize=8, fontname='Arial')
        axes[1].set_title("PR Curve", fontsize = 8, fontname = 'Arial')
        axes[1].set_xlim(-0.05, 1.05)
        axes[1].set_ylim(-0.05, 1.05)

        for label in axes[1].get_xticklabels() + axes[1].get_yticklabels():
            label.set_fontname('Arial')
            label.set_fontsize(8)

        pr_legend = axes[1].legend(loc="lower left", frameon=True)
        pr_legend.set_title(None)
        for text in pr_legend.get_texts():
            text.set_fontname('Arial')
            text.set_fontsize(8)

        # Final layout
        plt.tight_layout()

        # Save combined plot
        fig.savefig(os.path.join(self.output_dir, 'ROC_PR_curve_test.png'), dpi=600, bbox_inches="tight")

        # Close to free memory
        plt.close(fig)



    def histogram(self):
        """
        Plotting the histogram for predicted probabilities on test and unseen datasets.
        """
        # Global font settings
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 8

        # Create subplots: 1 row x 2 columns
        fig, axes = plt.subplots(1, 2, figsize=(8, 4), dpi=600)

        datasets = [('Test Set', self.pred_proba['VID']), ('Unseen Set', self.meta_unknown['infection_probability'])]

        for ax, (label, data) in zip(axes, datasets):
            sns.histplot(data, bins=30, color='darkblue', ax=ax)

            # Add threshold line if defined
            if self.threshold is not None:
                ax.axvline(x=self.threshold, color='red', linestyle='dashed', label=f'Threshold = {self.threshold:.2f}')
                ax.legend(loc='upper center', fontsize=8)

            # Styling
            ax.set_xlabel("Infection Probability", fontsize = 8)
            ax.set_ylabel("Count", fontsize = 8)
            ax.set_title(label, fontsize = 8)
            ax.grid(False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        plt.tight_layout()
        # Save the plot
        plt.savefig(os.path.join(self.output_dir, 'Infection_probability_histogram.png'), dpi=600, bbox_inches="tight")
        
        plt.close()

        
    def xgb_feature_importance(self):
        '''
        Visualize the feature importance of XGBClassifier in a grid layout.
        '''
        # Set font
        plt.rcParams['font.family'] = 'Arial'
        plt.rcParams['font.size'] = 8

        # Set up subplots: 1 row x 3 cols
        fig, axes = plt.subplots(3, 1, figsize=(3, 6), dpi=600)

        importance_types = ['weight', 'gain', 'cover']
        titles = {
            'weight': 'Frequency (Weight)',
            'gain': 'Average Gain',
            'cover': 'Average Coverage'
        }

        # Get booster from final estimator
        model = self.grid_result.best_estimator_.final_estimator_.get_booster()
        
        # Set your custom feature names
        model.feature_names = ['RF', 'SVM', 'KNN', 'GNB', 'LGR']
        
        # Plot feature importance for each type (weight, gain, cover)
        for ax, imp_type in zip(axes, importance_types):
            xgb.plot_importance(
                model,
                ax=ax,
                importance_type=imp_type,
                show_values=False,
                grid=False,
                color='darkblue'
            )

            # Format values to two decimal places by modifying the x-tick labels
            for tick in ax.get_xticklabels():
                tick.set_text(f'{float(tick.get_text()):.2f}')

            # Remove the title for the subplot
            ax.set_title('', fontsize=8)
            ax.set_xlabel(titles[imp_type], fontsize=8)
            ax.set_ylabel('Model', fontsize=8)  # Set y-axis label as 'Feature'
            
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)

        plt.tight_layout(pad=1.0)  # Reduce padding between subplots
        # Save the plot
        plt.savefig(os.path.join(self.output_dir, 'XGB_feature_importance.png'), dpi=600, bbox_inches="tight")
        plt.close()
    


        
    def evaluate_models_with_bootstrap(self, n_bootstrap=1000):
        """
        Evaluate multiple models using bootstrap sampling for classification metrics.

        Args:
            n_bootstrap: int, the number iteration for bootstrap.
        """
        # Define scoring functions
        metrics = {
            'Balanced Accuracy': balanced_accuracy_score,
            'Precision': precision_score,
            'Recall': recall_score,
            'Specificity': recall_score,
            'F1-score': f1_score,
            'ROC_AUC': roc_auc_score,
        }

        # Initialize structure to store bootstrapped scores
        bootstrap_scores = {
            metric: {model: [] for model in self.pred_proba}
            for metric in metrics
        }

        # Perform bootstrapping
        for _ in range(n_bootstrap):
            # Resample with replacement until both classes are present
            while True:
                indices = resample(np.arange(len(self.y_test)))
                y_true_boot = self.y_test[indices]
                if len(np.unique(y_true_boot)) == 2:
                    break

            for model, probs in self.pred_proba.items():
                y_prob = np.array(probs)[indices]
                y_pred = (y_prob >= 0.5).astype(int)

                for metric_name, metric_func in metrics.items():
                    if metric_name == 'ROC_AUC': 
                        score = metric_func(y_true_boot, y_prob)
                    elif metric_name == 'Balanced Accuracy':
                        score = metric_func(y_true_boot, y_pred)
                    elif metric_name == 'Specificity':
                        score = metric_func(y_true_boot, y_pred, pos_label = 0, zero_division = 0)
                    else:
                        score = metric_func(y_true_boot, y_pred, zero_division=0)
                    bootstrap_scores[metric_name][model].append(score)

        # Compute mean and CI
        results_dict = {}

        for metric_name, scores_by_model in bootstrap_scores.items():
            metric_df = pd.DataFrame(columns=['mean_score', 'ci_lower', 'ci_upper'])

            for model, scores in scores_by_model.items():
                mean_score = np.mean(scores)
                ci_lower = np.percentile(scores, 2.5)
                ci_upper = np.percentile(scores, 97.5)

                metric_df.loc[model] = [mean_score, ci_lower, ci_upper]

            results_dict[metric_name] = metric_df

        return results_dict


    def forest_plot(self):
        """
        Draw the forest plot. 
        """
        for key in self.ci_dict.keys():
            df = self.ci_dict[key]
            df = df.round(2)
            # Set global font
            matplotlib.rcParams['font.family'] = 'Arial'
            matplotlib.rcParams['font.size'] = 8
            
            # Create forest plot
            plot = EffectMeasurePlot(label=list(df.index), effect_measure=list(df['mean_score']), lcl=list(df['ci_lower']), ucl=list(df['ci_upper']))
            plot.labels(effectmeasure='Mean', center=max(df['mean_score']))
            plot.colors(pointshape="D", color='red')
            
            # Adjusted figure size
            fig_width, fig_height = 4, 1.5
            ci_range = df['ci_upper'].max() - df['ci_lower'].min()
            padding = 0.1 * ci_range
            max_value = round(df['ci_upper'].max() + 1.5 * padding, 2)
            min_value = round(max(df['ci_lower'].min() - padding, 0), 2)
            ax = plot.plot(figsize=(fig_width, fig_height), 
                           t_adjuster=0.06,
                           max_value=max_value, 
                           min_value=min_value,
                           text_size=8,
                           decimal=2,
                           size=1)
            
            # Title
            plt.title(f"{key}", loc="center", x=-0.82, y=1.045, fontsize=10, fontweight = 'bold')
            plt.suptitle("Model", x=0.08, y=0.913, fontsize=8)
            
            # Axis styling
            ax.set_xlabel("Lower Performance                        Higher Performance", fontsize=8)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            
            # Save the histogram of predicted probability
            plt.savefig(os.path.join(self.output_dir, f'{key}_forest.png'), dpi=600, bbox_inches="tight")  

            # Close the plot to free up memory
            plt.close()



    def combine_images_pillow(self):
        """
        Combine PNG images into a grid using PIL with added space between images.
        """
        # set the layout parameters
        grid_shape=(3, 2)
        padding=250
        
        # load image and get image attributes
        image_files = ['Balanced Accuracy_forest.png',
                       'Precision_forest.png',
                       'Recall_forest.png',
                       'Specificity_forest.png',
                       'F1-score_forest.png',
                       'ROC_AUC_forest.png']
        images = [Image.open(os.path.join(self.output_dir, f)) for f in image_files]
        rows, cols = grid_shape
        width, height = images[0].size

        # Create a blank canvas with extra space for padding between images
        combined_img = Image.new(
            'RGB', 
            ((cols * width) + (cols - 1) * padding, (rows * height) + (rows - 1) * padding), 
            (255, 255, 255)
        )

        # paste the image on the canvas
        for idx, img in enumerate(images):
            if idx >= rows * cols:
                break
            x_offset = (idx % cols) * (width + padding)
            y_offset = (idx // cols) * (height + padding)
            combined_img.paste(img, (x_offset, y_offset))

        # save the iamge
        combined_img.save(os.path.join(self.output_dir, 'Forest_plot_test.png'), dpi=(600, 600))
        
        # remove the separate images
        for image in image_files:
            os.remove(os.path.join(self.output_dir, image))
            
            
    def calibration_plot(self):
        """
        Draw the calibration file for all models
        """
        matplotlib.rcParams['font.family'] = 'Arial'
        matplotlib.rcParams['font.size'] = 8
        fig = plt.figure(figsize=(12, 10))
        gs = GridSpec(5, 2)
        colors = plt.get_cmap("Dark2")

        ax_calibration_curve = fig.add_subplot(gs[:2, :2])
        calibration_displays = {}
        for i, (clf, name) in enumerate(self.clf_list):
            if name == 'VID':
                X_test = self.X_test_norm[self.features]
            else:
                X_test = self.X_test_norm
            display = CalibrationDisplay.from_estimator(
                clf,
                X_test,
                self.y_test,
                n_bins=10,
                name=name,
                ax=ax_calibration_curve,
                color=colors(i),
            )
            calibration_displays[name] = display

        ax_calibration_curve.grid()
        ax_calibration_curve.set_title("Calibration Curves")

        # Add histogram
        grid_positions = [(2, 0), (2, 1), (3, 0), (3, 1), (4, 0), (4, 1)]
        for i, (_, name) in enumerate(self.clf_list):
            row, col = grid_positions[i]
            ax = fig.add_subplot(gs[row, col])

            ax.hist(
                calibration_displays[name].y_prob,
                range=(0, 1),
                bins=10,
                label=name,
                color=colors(i),
            )
            ax.set(title=name, xlabel="Mean predicted probability", ylabel="Count")

        plt.tight_layout()
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, f'Calibration_plot_test.png'), dpi=600, bbox_inches="tight")  

        # Close the plot to free up memory
        plt.close()