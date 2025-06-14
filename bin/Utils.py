import matplotlib.pyplot as plt
import matplotlib
from sklearn.metrics import confusion_matrix, roc_auc_score, ConfusionMatrixDisplay, roc_curve, RocCurveDisplay, PrecisionRecallDisplay, precision_score, recall_score, f1_score, balanced_accuracy_score
from boruta import BorutaPy
from args_dict import args_fs
from sklearn.utils import resample
import shap
import os
import xgboost as xgb
from sklearn.ensemble import RandomForestClassifier
import seaborn as sns
from zepid.graphics import EffectMeasurePlot
from matplotlib.gridspec import GridSpec
from sklearn.calibration import CalibratedClassifierCV, CalibrationDisplay
from PIL import Image
import numpy as np
import pandas as pd
from sklearn.utils import resample
from collections import Counter
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
        args_fs['estimator']['random_state'] = self.random_state
        args_fs['boruta']['random_state'] = self.random_state
        self.fs_estimator = RandomForestClassifier(**args_fs['estimator'])
        args_fs['boruta']['estimator'] = self.fs_estimator
        args_fs['boruta']['verbose'] = self.verbose
        args_fs['boruta']['max_iter'] = self.fs_iter
        boruta_selector = BorutaPy(**args_fs['boruta'])
        return boruta_selector

    def stratified_downsample(self, df, meta, target_size, random_state=42):
        """
        Stratified downsample based on class distribution.
        """
        grouped = []
        for label in meta['infection_status'].unique():
            # Filter both df and meta by label
            idx_label = meta[meta['infection_status'] == label].index
            df_label = df.loc[idx_label]
            
            # Determine number of samples for this label
            n_samples = int(target_size * (len(df_label) / len(meta)))
            
            # Resample (preserves index)
            df_sampled = resample(df_label, n_samples=n_samples, replace=False, random_state=random_state)
            
            # Sort sampled by sort_col
            df_sampled_sorted = df_sampled.loc[meta.loc[df_sampled.index].sort_values(by='infection_probability', ascending=False).index]
            
            grouped.append(df_sampled_sorted)
            
        return pd.concat(grouped)


    def get_stratified_background(self):
        """
        Create a stratified SHAP background sample, using double the size of the minority class.

        Parameters:
        - X: pd.DataFrame or np.ndarray, feature matrix
        - y: pd.Series or np.ndarray, target labels

        Returns:
        - background_X: pd.DataFrame, stratified background sample
        """
        df = self.X_train_norm[self.features].copy()
        df['__label__'] = self.y_train

        # Count instances per class
        class_counts = df['__label__'].value_counts()
        minority_class_count = class_counts.min()
        
        # Determine max_samples: 2 × minority class size
        if minority_class_count > 50:
            samples_per_class = 50
        else:
            samples_per_class = minority_class_count  

        stratified_samples = []

        for cls in class_counts.index:
            cls_samples = df[df['__label__'] == cls]
            sampled = resample(cls_samples,
                            replace=False,
                            n_samples=samples_per_class,
                            random_state=self.random_state)
            stratified_samples.append(sampled)

        background_df = pd.concat(stratified_samples).drop(columns='__label__')
        return background_df

    def shap_plot(self, X_target):
        """
        Draw the SHAP value beeswarm plot and heatmap for important genes in a 1x2 grid layout.
        """

        # Create SHAP explainer (for binary classification)
        background_df = self.get_stratified_background()
        masker = shap.maskers.Independent(background_df)
        explainer = shap.explainers.Permutation(
            self.grid_result.best_estimator_.predict_proba,
            masker,
            max_evals=2 * len(self.features) + 1
        )
        shap_values = explainer(X_target[self.features])

        # Take SHAP values for class 1 (malignant)
        shap_values_class1 = shap_values[..., 1]

        # Determine number of top features
        num_feature = min(20, len(self.features))

        # Get top N features by mean absolute SHAP
        mean_shap = np.abs(shap_values_class1.values).mean(axis=0)
        top_indices = np.argsort(mean_shap)[-num_feature:][::-1]
        shap_values = shap_values_class1[:, top_indices]
        X_test = X_target[self.features].iloc[:, top_indices]

        # Create a shared Explanation object
        explanation = shap.Explanation(
            values=shap_values.values,
            base_values=shap_values.base_values,
            data=X_test,
            feature_names=X_test.columns
        )

        # Create side-by-side plots (1 row, 2 columns)
        fig, axes = plt.subplots(1, 2, figsize=(18, 6))  # Adjust figsize if needed

        # Beeswarm plot on the left
        shap.plots.beeswarm(
            explanation,
            ax=axes[0],
            show=False,
            max_display=num_feature,
            plot_size=None
        )
        axes[0].set_xlabel('SHAP value', fontsize=12)
        axes[0].set_ylabel(f'Top {num_feature} genes', fontsize=12)
        axes[0].tick_params(axis='x', colors='black', labelsize=10)
        axes[0].tick_params(axis='y', colors='black', labelsize=10)

        # Heatmap on the right
        shap.plots.heatmap(
            explanation,
            ax=axes[1],
            show=False,
            max_display=num_feature,
            instance_order=np.arange(len(explanation))
        )
        
        # customized setting for heatmap
        axes[1].set_xlabel('Instances', fontsize=12)
        axes[1].set_ylabel(f'Top {num_feature} genes', fontsize=12)
        axes[1].tick_params(axis='x', colors='black', labelsize=10)
        axes[1].tick_params(axis='y', colors='black', labelsize=10)
        heatmap_cbar = fig.axes[-1]
        heatmap_cbar.set_ylabel("SHAP value", fontsize=12)
        
        # set the vertical line
        labels = [self.line_label[0], self.line_label[1], 100 - self.line_label[1], 100 - self.line_label[0]]
        self.line_label = []
        for label in labels:
            if label < 10:
                label = ' ' + str(label)
            else:
                label = str(label)
            self.line_label.append(label)
        for edge, label in zip(self.xticks[1:-1], 
                               [f'TP   PP\n{self.line_label[0]}%  {self.line_label[1]}%', 
                                '', 
                                f'PN   TN\n{self.line_label[2]}%  {self.line_label[3]}%']):
            axes[1].axvline(x=edge, color='gray', linestyle='--', linewidth=1)
            axes[1].text(
                edge,                     
                axes[1].get_ylim()[1],  
                label,
                ha='center',
                va='bottom',
                fontsize=10,
            )

        # Save combined plot
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, f"SHAP_plot_unseen.png"), dpi=600, bbox_inches='tight')
        plt.close()

    
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
        ax.set_xticklabels(disp.display_labels, fontsize=8, va='center')
        ax.set_yticklabels(disp.display_labels, fontsize=8, rotation = 90, va='center')

        # Keep the tick labels 
        ax.tick_params(axis='x', which='both', length=0, pad = 8)  # Removes x-axis ticks
        ax.tick_params(axis='y', which='both', length=0)  # Removes y-axis ticks
        
        # Set axis labels
        ax.set_xlabel('Prediction', fontsize=8, loc='center')
        ax.set_ylabel('Ground Truth', fontsize=8, loc='center')

        # Customize number annotations inside matrix
        for text in disp.text_.ravel():
            text.set_fontsize(8)

        # Colorbar font customization
        cbar = disp.im_.colorbar
        cbar.ax.tick_params(labelsize=8)
            
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

        axes[0].set_xlabel("False Positive Rate", fontsize=8)
        axes[0].set_ylabel("True Positive Rate", fontsize=8)
        axes[0].set_title("ROC Curve", fontsize = 8)
        axes[0].set_xlim(-0.05, 1.05)
        axes[0].set_ylim(-0.05, 1.05)

        for label in axes[0].get_xticklabels() + axes[0].get_yticklabels():
            label.set_fontsize(8)

        roc_legend = axes[0].get_legend()
        if roc_legend:
            roc_legend.set_title(None)
            for text in roc_legend.get_texts():
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

        axes[1].set_xlabel("Recall", fontsize=8)
        axes[1].set_ylabel("Precision", fontsize=8)
        axes[1].set_title("PR Curve", fontsize = 8)
        axes[1].set_xlim(-0.05, 1.05)
        axes[1].set_ylim(-0.05, 1.05)

        for label in axes[1].get_xticklabels() + axes[1].get_yticklabels():
            label.set_fontsize(8)

        pr_legend = axes[1].legend(loc="lower left", frameon=True)
        pr_legend.set_title(None)
        for text in pr_legend.get_texts():
            text.set_fontsize(8)

        # Final layout
        plt.tight_layout()

        # Save combined plot
        fig.savefig(os.path.join(self.output_dir, 'ROC_PR_curve_test.png'), dpi=600, bbox_inches="tight")

        # Close to free memory
        plt.close(fig)

        
    def xgb_feature_importance(self):
        '''
        Visualize the feature importance of XGBClassifier in a grid layout.
        '''
        # Set font
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
    


    def cv_bar_plot(self):
        '''
        Visualize the cross validation score with the bar plot.
        '''
        best_idx = self.grid_result.best_index_
        metrics = []
        means = []
        stds =  []
        for metric in self.scoring.keys():
            metrics.append(metric)
            means.append(self.grid_result.cv_results_[f'mean_test_{metric}'][best_idx])
            stds.append(self.grid_result.cv_results_[f'std_test_{metric}'][best_idx])
        
        # Compute dynamic ylim
        ymin = max(0, min([m - s for m, s in zip(means, stds)]) - 0.05)
        ymax = min(1.05, max([m + s for m, s in zip(means, stds)]) + 0.05)

        # Create figure
        plt.figure(figsize=(8, 5))
        bars = plt.bar(metrics, means, yerr=stds, capsize=8, color='skyblue', edgecolor='black')

        # Add text labels: mean (std) above error bar
        for bar, mean, std in zip(bars, means, stds):
            yval = bar.get_height() + std + 0.005  # position above error bar
            plt.text(
                bar.get_x() + bar.get_width() / 2,
                yval,
                f'{mean:.3f} (±{std:.3f})',
                ha='center',
                va='bottom',
                fontsize=10
            )

        # Titles and labels
        plt.ylim(ymin, ymax)
        plt.title('Cross-Validation Scores (Mean ± Std)', fontsize=12, fontweight='bold')
        plt.ylabel('Scores', fontsize=10, fontweight='bold')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'CV_bar.png'), dpi=600, bbox_inches="tight")
        plt.close()