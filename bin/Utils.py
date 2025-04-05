import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, roc_auc_score, ConfusionMatrixDisplay, roc_curve, RocCurveDisplay
from boruta import BorutaPy
from args_dict import args_fs
from sklearn.ensemble import RandomForestClassifier
import os
import xgboost
import seaborn as sns
import numpy as np
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
        Parameters:
        - pred: predictions of model on testing set
        - label: true class label of testing set

        Returns:
        - A confusion matrix plot with larger values.
        """
        cm = confusion_matrix(label, pred)

        # Create figure and axis
        fig, ax = plt.subplots(figsize=(4, 3))
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
        output_path = os.path.join(self.output_dir, 'Confusion_Matrix_test.png')
        plt.savefig(output_path, bbox_inches="tight", dpi=600)
        plt.close()

    
    def roc_plot(self, y_true, y_pred_proba):
        """
        Draw ROC plot.
        Args:
            y_true : np.array, the list of true label.
            y_pred_proba : np.array, the list of predicted probability.
        """
        # set the font type as 'Arial'
        plt.rcParams['font.family'] = 'Arial'

        fig, ax = plt.subplots(figsize=(4, 4))
        
        # Generate the ROC curve
        display = RocCurveDisplay.from_predictions(
            y_true,
            y_pred_proba,
            name="VID",
            color="blue",
            plot_chance_level=True,
        )

        ax = display.ax_
        

        # Customize the plot labels with Arial
        ax.set_xlabel("False Positive Rate", fontsize=9, fontname='Arial')
        ax.set_ylabel("True Positive Rate", fontsize=9, fontname='Arial')

        # Set custom limits to add spacing around the plot
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(-0.05, 1.05)
    
        # Set tick labels font
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontname('Arial')
            label.set_fontsize(8)
            
        # Customize legend
        legend = ax.get_legend()
        if legend:
            legend.set_title(None)  # Optional: remove legend title
            for text in legend.get_texts():
                text.set_fontname('Arial')
                text.set_fontsize(8)
        # Remove spines using seaborn
        #sns.despine(ax=ax)

        # Save the figure
        display.figure_.savefig(os.path.join(self.output_dir, 'ROC_Curve_test.png'), dpi=600, bbox_inches="tight")

        # Close the figure to free memory
        plt.close(display.figure_)

    def histogram(self, pred_proba, title = 'test'):
        """
        Plotting the histogram for predicted probability.
        Arguments:
            pred_proba: list or pd.DataFrame or np.array, the list of predicted probability
            title: str ('test' or 'unseen', default 'test'), the dataset to draw
        """
        plt.rcParams['font.family'] = 'Arial'
        # Plot the histogram
        plt.figure(figsize=(6, 4))
        sns.histplot(pred_proba, bins=30, color='black')
        
        # Add a vertical line at x = threshold
        if self.threshold:
            plt.axvline(x=self.threshold, color='red', linestyle='dashed', label=f'Threshold = {self.threshold}')
            plt.legend(loc='upper center')
        
        # Remove grid
        plt.grid(False)
        
        # remove the upd an right boundary
        ax = plt.gca()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Labels and title
        plt.xlabel("Infection Probability", fontweight = 'bold')
        plt.ylabel("Count", fontweight = 'bold')
        plt.title(f'Distribution of Infection Probability Predicted by {self.metamodel}', fontweight = 'bold')
        

        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, f'pred_proba_hist_{title}.png'), dpi=600, bbox_inches="tight")  
        
        # Close the plot to free up memory
        plt.close()
        
    def xgb_feature_importance(self):
        '''
        Visualize the feature importance of xgbclassifier.
        '''
        plt.rcParams['font.family'] = 'Arial'
        fig, ax = plt.subplots(figsize=(5, 4))

        xgboost.plot_importance(self.grid_result.best_estimator_, 
                            ax=ax,
                            grid=False,
                            importance_type='weight', 
                            title='Feature Importance of XGBoostClassifier',
                            show_values=True)

        # Remove the top and right edges
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'XGB_featrue_importance.png'), dpi=600, bbox_inches="tight")  

        # Close the plot to free up memory
        plt.close()

        
    def box_cv(self):
        plt.rcParams['font.family'] = 'Arial'
        df = self.val_cv_scores.iloc[[-1], :]

        # Step 1: Reshape the data
        # Filter only the split columns (e.g., test_accuracy and test_f1)
        split_columns = [col for col in df.columns if col.startswith('split')]
        reshaped_data = df[split_columns].melt(var_name="Metric_Fold", value_name="Value")
        
        # Extract metric names from column names
        reshaped_data['Metric'] = reshaped_data['Metric_Fold'].apply(lambda x: "_".join(x.split('_')[2:])) 
        reshaped_data['Fold'] = reshaped_data['Metric_Fold'].apply(lambda x: x.split('_')[0])  
        
        # Step 2: Draw the boxplot
        plt.figure(figsize=(5, 4))
        sns.boxplot(data=reshaped_data, x='Metric', y='Value', palette="Set2")
        
        # Step 3: Overlay data points
        sns.stripplot(data=reshaped_data, x='Metric', y='Value', color='black', size=4, jitter=True, alpha=0.6)
        
        # Customizing the plot
    #    plt.title("Boxplot of cross validation scores", fontsize=14)
        plt.xlabel("Metrics", fontsize=8)
        plt.ylabel("Values", fontsize=8)
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.tick_params(axis='both', which='major', labelsize=8)
        # Remove the right and top spines
        sns.despine()
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'cv_box.png'), dpi=600, bbox_inches="tight")  

        # Close the plot to free up memory
        plt.close()