import shutil

import numpy as np
import pandas as pd
import os
import time
from sklearn.metrics import auc
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import RocCurveDisplay
from sklearn.feature_selection import SelectFromModel
from matplotlib import pyplot as plt

"""
1. This file is used to predict whether the amino acid sequence file output by getPossibleIFNSeqs.py corresponds to IFN. For detailed procedures, refer to the article. 

2. Parameters to be modified:

    'neg_sample_path': Directory for the negative samples (neighborhood-genes for IFN) used to train the model; 
     these files are compiled by `collectPosNegSamples.py`.
    
    'pos_sample_path': Directory for the positive samples (IFN) used to train the model.
    
    'root_path': Root directory; primarily used for outputting results and locating the iFeature folder (implementing the APAAC method).
    
    'query_path': Directory for the samples to be predicted (amino acid sequence files output by getPossibleIFNSeqs.py.py).
"""

neg_sample_path = 'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\IFNSCOPE/neg_samples.fasta'
pos_sample_path = 'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\IFNSCOPE/pos_samples.fasta'
root_path = 'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs'
query_path = 'E:\Study\Doctor\Projects\IFN\代码\mfb\IFNSeqs\queryByClustersOutput'



classifier = LogisticRegression(solver='sag', penalty='l2', max_iter=50000, random_state=919)


def mkdir(path, folder_name):
    folder_path = f'{path}/{folder_name}'
    exist = os.path.exists(folder_path)
    if not exist:  # Check if the directory exists; if it does not, create it as a directory.
        os.makedirs(folder_path)
    return folder_path


def get_filenames(path):
    names = os.listdir(path)
    return names


def collect_all_samples(pos_path, neg_path, rpath):
    with open(pos_path) as f:
        pos_samples = f.read()
    pos_num = pos_samples.count('>')
    with open(neg_path) as f:
        neg_samples = f.read()
    neg_num = neg_samples.count('>')
    result_path = mkdir(rpath, 'IFNSCOPE')
    all_samples_path = f'{result_path}/all_samples.txt'
    with open(all_samples_path, 'w') as f:
        f.write(pos_samples + neg_samples)
    return result_path, pos_num, neg_num


def feature_extraction(result_path, qpath, rpath):
    ifeature_path = f'{rpath}/iFeature'
    all_samples_path = f'{result_path}/all_samples.txt'
    os.system(f'python {ifeature_path}/iFeature.py --file {all_samples_path} --type APAAC ')
    os.rename('encoding.tsv', 'train_APAAC.txt')
    shutil.move('train_APAAC.txt', f'{result_path}')
    time.sleep(2)
    possible_IFNs_path = f'{qpath}/possible_IFNs_summary.fasta'
    os.system(f'python {ifeature_path}/iFeature.py --file {possible_IFNs_path} --type APAAC')
    os.rename('encoding.tsv', 'predict_APAAC.txt')
    shutil.move('predict_APAAC.txt', f'{result_path}')

    train_path = f'{result_path}/train_APAAC.txt'
    predict_path = f'{result_path}/predict_APAAC.txt'
    return train_path, predict_path


def get_x_and_y(train_path, predict_path, pos_num, neg_num):
    x_df = pd.read_csv(train_path, sep='\t')
    x_pred_df = pd.read_csv(predict_path, sep='\t')
    x = x_df.iloc[:, 1:]
    features = x.columns
    x_pred = x_pred_df.iloc[:, 1:]
    sc = StandardScaler()
    x = sc.fit_transform(x)
    x_pred = sc.transform(x_pred)
    y = []
    for n in range(pos_num):
        y.append(1)
    for n in range(neg_num):
        y.append(0)
    y = np.array(y)
    return x, x_pred, y, features


def feature_selection(x, x_predicted, y, features):
    # This method is based on the SelectFromModel method from sklearn.feature_selection.
    # # It is model-based and filters features with correlation coefficients below the threshold.

    # Before dimensionality reduction
    mean_scroe = cross_val_score(classifier, x, y, cv=5, scoring='f1').mean()
    print(f'mean_f1 before feature selection: {mean_scroe}')
    # After dimensionality reduction
    selector = SelectFromModel(estimator=classifier).fit(x, y)
    x_selected = selector.transform(x)
    x_pred_selected = selector.transform(x_predicted)
    coef = selector.estimator_.coef_[0]
    isSelected = selector.get_support()
    plot_feature_importance(features, coef, isSelected)
    mean_scroe2 = cross_val_score(classifier, x_selected, y, cv=5, scoring='f1').mean()
    print(f'mean_f1 after feature selection: {mean_scroe2}')
    return x_selected, x_pred_selected


def plot_feature_importance(features, coef, isSelected):
    features_omit = []
    for feature in features:
        feature_omit = feature.split('.', 1)[1]
        if len(feature_omit) > 2:
            split_list = feature_omit.split('.', 1)
            if 'phobi' in split_list[0]:
                feature_omit = f'phobi.{split_list[1]}'  # split_list[1] is the distance of aa
            else:
                feature_omit = f'phili.{split_list[1]}'
        features_omit.append(feature_omit)
    array = np.array(features_omit)

    importances = pd.DataFrame(data={
        'Attribute': array,
        'Importance': coef
    })
    importances = importances.sort_values(by='Importance', ascending=False)

    color_list = []
    for idx in importances.index:
        if isSelected[idx]:
            # Selected features
            color_list.append('#FFC0CB')
        else:
            color_list.append('#D3D3D3')

    plt.figure(figsize=(12, 6), dpi=100)
    plt.bar(x=importances['Attribute'], height=importances['Importance'], color=color_list)
    # plt.title('Feature dimensional reduction result', size=18)
    # plt.xlabel('Features', fontsize=18)
    plt.xticks(rotation=60, rotation_mode='anchor', ha='right', fontsize=8.5)
    # plt.ylabel('Coefficients', fontsize=18)
    plt.yticks(fontsize=12)
    # plt.show()
    # plt.savefig("Dimensional reduction.svg", dpi=400, format="svg")


def model_evaluation(x, y):
    mean_f1 = cross_val_score(classifier, x, y, cv=5, scoring='f1').mean()
    mean_accuracy = cross_val_score(classifier, x, y, cv=5, scoring='accuracy').mean()
    mean_auc = draw_roc_plot(classifier, x, y)
    print(f'mean_f1: {mean_f1}')
    print(f'mean_accuracy: {mean_accuracy}')
    print(f'mean_auc: {mean_auc}')


def draw_roc_plot(classifier, x, y):
    cv = StratifiedKFold(n_splits=5)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    i = 0
    fig, ax = plt.subplots()
    for i, (train, test) in enumerate(cv.split(x, y)):
        classifier.fit(x[train], y[train])
        viz = RocCurveDisplay.from_estimator(classifier, x[test], y[test])
        interp_tpr = np.interp(mean_fpr, viz.fpr, viz.tpr)
        interp_tpr[0] = 0.0
        tprs.append(interp_tpr)
        aucs.append(viz.roc_auc)
    ax.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey',
            label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    ax.plot(mean_fpr, mean_tpr, color='red',
            label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
            lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    ax.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                    label=r'$\pm$ 1 std. dev.')

    ax.set(xlim=[-0.05, 1.05], ylim=[-0.05, 1.05],
           title="Receiver operating characteristic curve")
    plt.xlabel('False positive rate', fontsize=20)
    plt.ylabel('True positive rate', fontsize=20)
    ax.legend(loc="lower right")
    # plt.show()
    return mean_auc


def train_and_predict(x, x_pred, y):
    classifier.fit(x, y)
    y_pred = classifier.predict(x_pred)
    pos_num = 0
    neg_num = 0
    for result in y_pred:
        if result:
            pos_num += 1
        else:
            neg_num += 1
    print(f'   Prediction complete, pos_num: {pos_num}, neg_num: {neg_num}')
    return y_pred


def plot_decision_boundary(classifier, x, y):
    plt.show()


def write_by_model_result(qpath, y):
    IFN_contents = ''
    possible_IFNs_path = f'{qpath}/possible_IFNs_summary.fasta'
    with open(possible_IFNs_path) as f:
        f_lines = f.readlines()
    idx = 0
    for n, line in enumerate(f_lines):
        if line.startswith('>'):
            if y[idx]:
                IFN_contents += line
                IFN_contents += f_lines[n + 1]
                IFN_contents += '\n'
            idx += 1
    predicted_IFNs_path = f'{qpath}/predicted_IFNs_summary.fasta'
    with open(predicted_IFNs_path, 'w') as f:
        f.write(IFN_contents)


if __name__ == "__main__":
    model_result_path, pos_sample_num, neg_sample_num = collect_all_samples(pos_sample_path, neg_sample_path, root_path)
    train_sample_path, predict_sample_path = feature_extraction(model_result_path, query_path, root_path)
    x, x_predicted, y, features = get_x_and_y(train_sample_path, predict_sample_path, pos_sample_num, neg_sample_num)
    x_selected, x_pred_selected = feature_selection(x, x_predicted, y, features)
    model_evaluation(x_selected, y)
    y_predicted = train_and_predict(x_selected, x_pred_selected, y)
    write_by_model_result(query_path, y_predicted)
