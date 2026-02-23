import os
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
from catboost import CatBoostClassifier
import lightgbm as lgb
import matplotlib.pyplot as plt
import shap

from feature_utils import extract_features

FEATURE_NAMES = [
    "gc_2_6",
    "gc_9_12",
    "gc_13_17",
    "gc_9_17",
    "tm_2_6_plus_9_17",
    "tt_count_2_6_plus_9_17",
    "aa_count_2_6_plus_9_17",
    "actt_count_2_6",
    "ctt_count_9_12",
    "tc_count_27_31",
    "is_pos35_c",
    "is_pos36_c",
    "ttc_count_22_26",
    "tc_count_22_26",
    "ttt_count_13_17",
    "tt_count_13_17",
    "is_pos6_gc",
    "is_pos9_gc",
    "is_pos17_gc",
    "is_pos22_gc",
    "is_pos31_gc",
    "is_pos35_gc",
]

def get_sequence_array(df):
    candidate_cols = ['Seqs', 'Seq', 'sg_sequence', 'Sequence']
    for col in candidate_cols:
        if col in df.columns:
            return df[col].astype(str).values

    if df.shape[1] >= 1:
        return df.iloc[:, -1].astype(str).values

    raise KeyError(
        f"No sequence column found. Available columns: {list(df.columns)}. "
        "Expected one of ['Seqs', 'Seq', 'sg_sequence', 'Sequence'] or at least 1 column for fallback to the last column."
    )


def main():

    file_path = "train_data/Screening_variants_by_lfc.xlsx"
    low_df = pd.read_excel(file_path, sheet_name='LowActivation')
    high_df = pd.read_excel(file_path, sheet_name='HighActivation')

    low_seqs = get_sequence_array(low_df)
    high_seqs = get_sequence_array(high_df)

    all_seqs = np.concatenate((low_seqs, high_seqs), axis=0)
    labels = []
    for i in range(low_seqs.size):
        labels.append(0)
    for i in range(high_seqs.size):
        labels.append(1)

    # Extract features
    features = []
    for seq in all_seqs:
        features.append(extract_features(seq))
    X = np.array(features)
    X_df = pd.DataFrame(X, columns=FEATURE_NAMES)

    Y = labels

    # Split data
    X_train, X_test, y_train, y_test = train_test_split(X_df, Y, test_size=0.2, stratify=Y)

    print("=" * 60)
    print("Boosting Models Comparison")
    print("=" * 60)

    results = {}

    # 1. XGBoost
    print("\n1. Training XGBoost...")
    xgb_model = XGBClassifier(max_depth=2, n_estimators=100)
    xgb_model.fit(X_train, y_train)

    # 2. CatBoost
    print("\n2. Training CatBoost...")
    cat_model = CatBoostClassifier(max_depth=2, iterations=100)
    cat_model.fit(X_train, y_train)

    # 3. LightGBM
    print("\n3. Training LightGBM...")
    lgb_model = lgb.LGBMClassifier(max_depth=2, n_estimators=100)
    lgb_model.fit(X_train, y_train)

    # Summary
    # Save comparison chart
    output_dir = "results/boosting_outputs"
    os.makedirs(output_dir, exist_ok=True)

    # Save feature importances for each model
    print("\n" + "=" * 60)
    print("Feature Importance (Top 10)")
    print("=" * 60)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

    models_dict = {
        'XGBoost': xgb_model,
        'CatBoost': cat_model,
        'LightGBM': lgb_model
    }

    for idx, (model_name, model) in enumerate(models_dict.items()):
        if model_name == 'XGBoost':
            importances = xgb_model.feature_importances_
        elif model_name == 'CatBoost':
            importances = cat_model.get_feature_importance()
        else:  # LightGBM
            importances = lgb_model.feature_importances_

        importance_df = pd.DataFrame({
            'Feature': FEATURE_NAMES,
            'Importance': importances
        }).sort_values('Importance', ascending=False).head(10)

        axes[idx].barh(importance_df['Feature'], importance_df['Importance'], color=colors[idx], alpha=0.7)
        axes[idx].set_xlabel('Importance')
        axes[idx].set_title(f'{model_name} - Top 10 Features')
        axes[idx].invert_yaxis()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "feature_importance_comparison.png"), dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Feature importance chart saved to: {os.path.join(output_dir, 'feature_importance_comparison.png')}")

    # SHAP Beeswarm plots
    print("\n" + "=" * 60)
    print("Generating SHAP Beeswarm Plots")
    print("=" * 60)

    shap_models = {
        'XGBoost': xgb_model,
        'CatBoost': cat_model,
        'LightGBM': lgb_model
    }

    for model_name, model in shap_models.items():
        print(f"\nGenerating SHAP plot for {model_name}...")
        try:
            explainer = shap.TreeExplainer(model)
            shap_values = explainer.shap_values(X_train)

            if isinstance(shap_values, list):
                shap_values = shap_values[1] if len(shap_values) > 1 else shap_values[0]

            plt.figure(figsize=(12, 8))
            shap.summary_plot(shap_values, X_train, show=False)
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir, f"shap_beeswarm_{model_name.lower()}.png"), dpi=300, bbox_inches="tight")
            plt.close()
            print(f"   SHAP beeswarm plot saved: shap_beeswarm_{model_name.lower()}.png")
        except Exception as e:
            print(f"   Warning: Could not generate SHAP plot for {model_name}: {str(e)}")

    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()
