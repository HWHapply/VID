import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, roc_auc_score, ConfusionMatrixDisplay, roc_curve
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
        
        fig, ax = plt.subplots(figsize=(10, 6))
        disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=['uninfected', 'infected'])
        
        disp.plot(include_values=True, cmap='Blues', values_format='d', ax=ax)
        
        # Increase the font size of the display labels
        ax.set_xticklabels(disp.display_labels, fontsize=14)
        ax.set_yticklabels(disp.display_labels, fontsize=14)
        
        # Increase the font size of the annotations
        for text in disp.text_.ravel():
            text.set_fontsize(20)
        
        plt.ylabel('Ground Truth', fontsize=16)
        plt.xlabel('Prediction', fontsize=16)
        plt.title(f'Confusion Matrix', fontsize=20)
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, 'Confusion_Matrix_test.png'))  
        
        # Close the plot to free up memory
        plt.close()

    def roc_plot(self, y_true, y_pred_prob):
        """
        Parameters:
        y_true: array-like, shape (n_samples,)
            True binary labels.
        y_pred_prob: array-like, shape (n_samples,)
            Target scores, can either be probability estimates of the positive class, confidence values, or non-thresholded measure of decisions.
        """

        # Compute ROC curve and ROC area
        fpr, tpr, _ = roc_curve(y_true, y_pred_prob)
        roc_auc = roc_auc_score(y_true, y_pred_prob)
        
        # Plotting
        plt.figure(figsize=(6, 6))
        plt.plot(fpr, tpr, color="b", label=f'ROC curve (AUC = {roc_auc:.2f})', lw=2, alpha=0.8)
        plt.plot([0, 1], [0, 1], linestyle="--", lw=2, color="r", label="Chance", alpha=0.8)
        
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xlabel('False Positive Rate', fontsize = 10)
        plt.ylabel('True Positive Rate', fontsize = 10)
        plt.title(f'ROC curve of {self.metamodel} on testing set', fontsize = 10)
        plt.legend(loc="lower right")

        # Save the plot to a directory
        plt.savefig(os.path.join(self.output_dir, 'ROC_Curve_test.png'))  # Replace with your desired directory and filename
        
        # Close the plot to free up memory
        plt.close()

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
        
        # Save the histogram of predicted probability
        plt.savefig(os.path.join(self.output_dir, f'pred_proba_hist_{title}.png'))  
        
        # Close the plot to free up memory
        plt.close()
        
    def xgb_feature_importance(self):
        '''
        Visualize the feature importance of xgbclassifier.
        '''
        fig, ax = plt.subplots(figsize=(4, 3))

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