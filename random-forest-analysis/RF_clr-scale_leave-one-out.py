import pandas as pd
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, roc_curve, auc, confusion_matrix, accuracy_score
from sklearn.impute import SimpleImputer
from sklearn.feature_selection import VarianceThreshold
from sklearn.preprocessing import StandardScaler, LabelEncoder
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

# Argument parser for input parameters
parser = argparse.ArgumentParser(description="Random Forest with Leave-One-Out Cross-Validation")
parser.add_argument("--threads", type=int, default=1, help="Number of threads for parallel processing")
parser.add_argument("--k_features", type=int, default=50, help="Number of top features to select")
parser.add_argument("--input_file", type=str, required=True, help="Input data file path")
args = parser.parse_args()

# Load data
data = pd.read_csv(args.input_file, sep='\t')

# Separate features and labels
X = data.drop(columns=['Sample', 'Group','Response', 'diarrhea_response', 'stomachache_response'])
y = data['stomachache_response']


# Impute missing values
imputer = SimpleImputer(strategy='median')
X_imputed = imputer.fit_transform(X)

# Filter low-variance features
selector = VarianceThreshold(threshold=0.01)
X_filtered = selector.fit_transform(X_imputed)
filtered_support = selector.get_support()
X_filtered_columns = np.array(X.columns)[filtered_support]

# Save filtered features with names
X_filtered_df = pd.DataFrame(X_filtered, columns=X_filtered_columns)
X_filtered_df.to_csv("filtered_features_with_names.tsv", sep="\t", index=False)

# Encode labels
label_encoder = LabelEncoder()
y_encoded = label_encoder.fit_transform(y)


# CLR Transformation
def clr_transform(matrix, pseudo_count=1e-6):
    matrix += pseudo_count
    geometric_mean = np.exp(np.mean(np.log(matrix), axis=1))
    clr_matrix = np.log(matrix / geometric_mean[:, np.newaxis])
    return clr_matrix


X_scaled = clr_transform(X_filtered)

# Calculate initial feature importance to select top k_features
model_rf = RandomForestClassifier(
    n_estimators=500,
    random_state=42,
    class_weight="balanced",
    n_jobs=args.threads
)
model_rf.fit(X_scaled, y_encoded)
feature_importances = model_rf.feature_importances_
top_k_indices = np.argsort(feature_importances)[::-1][:args.k_features]
top_features = X_filtered_columns[top_k_indices]
top_importances = feature_importances[top_k_indices]

# Save top features and their importances to a CSV file
feature_importance_df = pd.DataFrame({"Feature": top_features, "Importance": top_importances})
feature_importance_df.to_csv("top_feature_importances_loocv.csv", index=False)
print(f"Top {args.k_features} features and their importances saved to 'top_feature_importances_loocv.csv'")

# Reduce X_scaled to top k_features
X_top_k = X_scaled[:, top_k_indices]

# Leave-One-Out Cross-Validation
loo = LeaveOneOut()
y_true = []
y_pred = []

print(f"Performing Leave-One-Out Cross-Validation with {len(X_top_k)} samples using top {args.k_features} features...")

for train_index, test_index in loo.split(X_top_k):
    X_train, X_test = X_top_k[train_index], X_top_k[test_index]
    y_train, y_test = y_encoded[train_index], y_encoded[test_index]

    model_rf.fit(X_train, y_train)
    y_true.append(y_test[0])
    y_pred.append(model_rf.predict(X_test)[0])

# Calculate accuracy
accuracy = accuracy_score(y_true, y_pred)
print(f"Leave-One-Out Cross-Validation Accuracy using top {args.k_features} features: {accuracy:.4f}")

# Confusion Matrix
conf_mat = confusion_matrix(y_true, y_pred)
plt.figure(figsize=(8, 6))
sns.heatmap(conf_mat, annot=True, fmt='d', cmap='Blues', xticklabels=np.unique(y), yticklabels=np.unique(y))
plt.title('Confusion Matrix (LOOCV)')
plt.xlabel('Predicted Labels')
plt.ylabel('True Labels')
plt.savefig("confusion_matrix_loocv.pdf")
plt.close()

# ROC Curve for binary classification
y_pred_proba = np.zeros(len(y_encoded))

for i, (train_index, test_index) in enumerate(loo.split(X_top_k)):
    X_train, X_test = X_top_k[train_index], X_top_k[test_index]
    y_train, y_test = y_encoded[train_index], y_encoded[test_index]

    model_rf.fit(X_train, y_train)
    y_pred_proba[test_index] = model_rf.predict_proba(X_test)[:, 1]

fpr, tpr, thresholds = roc_curve(y_encoded, y_pred_proba)
roc_auc = auc(fpr, tpr)

# Plot ROC curve
plt.figure(figsize=(12, 9))
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f"ROC curve (AUC = {roc_auc:.2f})")
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label="Random Guess")
plt.title("ROC Curve with LOOCV")
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.legend(loc="lower right")
plt.savefig("binary_roc_curve_loocv.pdf")
plt.close()

print("LOOCV with top features and visualizations complete.")
