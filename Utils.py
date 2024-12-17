import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, roc_auc_score, ConfusionMatrixDisplay, roc_curve, RocCurveDisplay
from boruta import BorutaPy
from args_dict import args_fs
from sklearn.ensemble import RandomForestClassifier
import os
import xgboost
import umap
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
        
        fig, ax = plt.subplots(figsize=(8, 6))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['uninfected', 'infected'])
        
        disp.plot(include_values=True, cmap='Blues', values_format='d', ax=ax)
        
        # Increase the font size of the display labels
        ax.set_xticklabels(disp.display_labels, fontsize=14)
        ax.set_yticklabels(disp.display_labels, fontsize=14)
        
        # Increase the font size of the annotations
        for text in disp.text_.ravel():
            text.set_fontsize(20)
        
        # set label and title format
        plt.ylabel('Ground Truth', fontsize=16)
        plt.xlabel('Prediction', fontsize=16)
        plt.title(f'Confusion Matrix', fontsize=20)
        
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'Confusion_Matrix_test.png'), bbox_inches='tight', dpi=300)  
        
        # Close the plot to free up memory
        plt.close()

    
    def roc_plot(self, y_true, y_pred_proba):
        # Generate the ROC curve
        display = RocCurveDisplay.from_predictions(
            y_true,
            y_pred_proba,
            name=f"positive vs negative",
            color="darkorange",
            plot_chance_level=True,
        )
        # Customize the plot labels and title
        _ = display.ax_.set(
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title=f"ROC curve",
        )
        # Remove spines using seaborn
        sns.despine(ax=display.ax_)
        display.figure_.savefig(os.path.join(self.output_dir, 'ROC_Curve_test.png'), dpi=300, bbox_inches="tight")

        # Close the figure to free memory
        plt.close(display.figure_)

    def histogram(self, pred_proba, title = 'test'):
        """
        Plotting the histogram for predicted probability.
        Arguments:
            pred_proba: list or pd.DataFrame or np.array, the list of predicted probability
            title: str ('test' or 'unseen', default 'test'), the dataset to draw
        """
        plt.figure()
        plt.hist(pred_proba, bins=10, range=(0, 1), color='skyblue', edgecolor='black')

        # Adding labels and title
        plt.xlabel('Probability')
        plt.ylabel('Frequency')
        plt.title(f'Distribution of Predicted Probability of {self.metamodel}')
        
        # remove the top and right edges
        sns.despine()
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, f'pred_proba_hist_{title}.png'))  
        
        # Close the plot to free up memory
        plt.close()
        
    def xgb_feature_importance(self):
        '''
        Visualize the feature importance of xgbclassifier.
        '''
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
        plt.savefig(os.path.join(self.output_dir, 'XGB_featrue_importance.png'))  

        # Close the plot to free up memory
        plt.close()


    def umap_plot(self):
        # Initialize UMAP model
        reducer = umap.UMAP()

        # Fit and transform the data
        embedding = reducer.fit_transform(self.data_df_processed)

        # Create a scatter plot of the UMAP embedding
        plt.figure(figsize=(8, 8))
        sns.scatterplot(
            x=embedding[:, 0], y=embedding[:, 1],
            hue=self.meta_df['infection_status'],
            palette=sns.color_palette("hsv", len(set(self.meta_df['infection_status']))),
            legend="full",
            alpha=0.7
        )
        plt.title('UMAP projection')
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'Umap_infection_status.png'))  

        # Close the plot to free up memory
        plt.close()
        
    def box_cv(self):
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
        plt.xlabel("Metrics", fontsize=12)
        plt.ylabel("Values", fontsize=12)
        plt.xticks(rotation=45)
        plt.tight_layout()
        # Remove the right and top spines
        sns.despine()
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'cv_box.png'))  

        # Close the plot to free up memory
        plt.close()