## this tool is the core function of prediction analysis
## author: taozhou
## email: zhou.tao@genecast.com.cn

import matplotlib as mpl
mpl.use('Agg')
import pandas as pd
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from scipy.stats import pearsonr
from sklearn.svm import LinearSVC, SVC
from sklearn.linear_model import LassoLarsIC, LogisticRegression, RandomizedLasso
from sklearn.feature_selection import RFE
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
import itertools
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")


class MethodException(Exception):
    pass


def loop_choose_svm(data, group, n):
    ## 循环遍历所有组合 进行预测 接受测试样本数量n
    a = group.columns[0]
    b = group.columns[-1]
    accuracy = defaultdict(list)
    for a_test in itertools.combinations(group[a], n):
        for b_test in itertools.combinations(group[b], n):
            a_train = [i for i in group[a] if i not in a_test]
            b_train = [i for i in group[b] if i not in b_test]
            feature_gene = feature_select(data, list(a_train), list(b_train))
            X, y = get_train_test_data(data, list(a_train), list(b_train), feature_gene)
            X1, y1 = get_train_test_data(data, list(a_test), list(b_test), feature_gene)
            for sample, re in zip(list(a_test) + list(b_test), svm(X, y, X1, y1)):
                accuracy[sample].append(re)
    accuracy = {key: value.count(True)/len(value) for key, value in accuracy.items()}
    return accuracy


def get_train_test_data(data, group1, group2, feature_gene):
    ## 获取训练及测试数据
    group1_data = data.ix[feature_gene][group1]
    group2_data = data.ix[feature_gene][group2]
    X = pd.concat([group1_data, group2_data], axis=1).T.as_matrix()
    y = [-1] * len(group1) + [1] * len(group2)
    return X, y


def feature_select(data, group1, group2, pval=0.05, method="logistic", criterion='aic', penalty="l2", \
                   C=1.0, threshold=0):
    ## 选择特征基因，支持多种选择方法， 默认采用logistic回归拟合权重系数
    regression_function = ["Lasso", "logistic", "RandomizedLasso"]
    if method not in regression_function:
        pearsonr_target = [-1] * len(group1) + [1] * len(group2)
        result = {"gene": data.index, "a_vs_b": []}
        for x, y in zip(data[group1].as_matrix(), data[group2].as_matrix()):
            if method == "wilcox":
                z, p = ranksums(x, y)
            elif method == "pearsonr":
                cor, p = pearsonr(pearsonr_target, list(x) + list(y))
            else:
                raise MethodException('must offer a right feature selection method')
            result["a_vs_b"].append(p)
        result = pd.DataFrame(result)
        result = result.loc[result["a_vs_b"] <= pval]["gene"]
        if len(result) == 0:
            raise MethodException("the number of feature gene is 0, please choose your fsm or other parameter")
        return list(result)
    else:
        sc = StandardScaler()
        X, y = get_train_test_data(data, group1, group2, data.index)
        sc.fit(X)
        X_train_std = sc.transform(X)
        genes = data.index
        if method == "Lasso":
            clf = LassoLarsIC(criterion=criterion)
        elif method == "logistic":
            clf = LogisticRegression(penalty=penalty, C=C)
        elif method == "RandomizedLasso":
            clf = RandomizedLasso()
        else:
            raise MethodException('must offer a right feature selection method')
        clf.fit(X_train_std, y)
        feature_gene = [gene for gene, value in zip(genes, *clf.coef_) if abs(value) >= threshold]
        if len(feature_gene) == 0:
            raise MethodException("the number of feature gene is 0, please choose your fsm or other parameter")
        return feature_gene


def evaluate_model(data, group1, group2, feature_gene, name=None, method="LinearSVC", C=1, n_folds=5):
    X, y = get_train_test_data(data, group1, group2, feature_gene)
    if method == "LinearSVC":
        model = LinearSVC()
    elif method == "SVC":
        model = SVC(kernel='linear', probability=True)
    elif method == "logistic":
        model = LogisticRegression(C=C, random_state=0)
    else:
        raise MethodException('must offer a right prediction method')
    plot_ROC(X, y, classifier=model, n_folds=n_folds, name=name)


def prediction(X, y, X1, y1, method="LinearSVC", C=1, n_folds=5):
    ## 支持向量机，支持多种子分类器，默认采用线性支持向量分类
    if method == "LinearSVC":
        model = LinearSVC()
    elif method == "SVC":
        model = SVC(kernel='linear', probability=True)
    elif method == "logistic":
        model = LogisticRegression(C=C, random_state=0)
    else:
        raise MethodException('must offer a right prediction method')
    model.fit(X, y)
    expected = y1
    predicted = model.predict(X1)
    plot_ROC(X, y, classifier=model, n_folds=n_folds, method=method)
    result = [True if i == k else False for i, k in zip(expected, predicted)]
    return result


def plot_ROC(X, y, classifier=None, n_folds=10, name=None, method=None):
    cv = StratifiedKFold(y, n_folds=n_folds)
    classifier = classifier
    mean_tpr = 0.0
    mean_fpr = np.linspace(0, 1, len(y))
    y = np.array(y)
    for i, (train, test) in enumerate(cv):
        probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        mean_tpr += interp(mean_fpr, fpr, tpr)
        mean_tpr[0] = 0.0
    mean_tpr /= len(cv)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.plot(mean_fpr, mean_tpr, color='darkorange', lw=2, label='Mean ROC (area = %0.2f)' % mean_auc)
    plt.xlim([-0.00, 1.00])
    plt.ylim([-0.00, 1.00])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title("ROC curve of %s (AUC = %.3f)" % (method, mean_auc))
    plt.savefig(name + "_ROC.png")
    plt.close()
    

if __name__ == "__main__":
    ## 主程序 返回并打印每个样本的准确率
    group = pd.read_table("group.txt")
    panel_gene = pd.read_table("sample_median_log2.txt", index_col=0)
    accuracy = loop_choose_svm(panel_gene, group, 1)
    result = pd.DataFrame(list(accuracy.values()), index=accuracy.keys(), columns=["accuracy"])
    print(result)

