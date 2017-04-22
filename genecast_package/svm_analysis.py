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
from scipy.stats import fisher_exact
from sklearn.svm import LinearSVC, SVC
from sklearn.linear_model import LassoLarsIC, LogisticRegression, RandomizedLasso, Lasso
from sklearn.feature_selection import RFE
from sklearn.feature_selection import VarianceThreshold
from sklearn.ensemble import RandomForestRegressor
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc
import itertools
from collections import defaultdict
import warnings
warnings.filterwarnings("ignore")


class MethodException(Exception):
    pass


class LR(LogisticRegression):
    def __init__(self, threshold=0.01, dual=False, tol=1e-4, C=1.0,
                 fit_intercept=True, intercept_scaling=1, class_weight=None,
                 random_state=None, solver='liblinear', max_iter=100,
                 multi_class='ovr', verbose=0, warm_start=False, n_jobs=1):

        # 权值相近的阈值
        self.threshold = threshold
        LogisticRegression.__init__(self, penalty='l1', dual=dual, tol=tol, C=C,
                                    fit_intercept=fit_intercept, intercept_scaling=intercept_scaling,
                                    class_weight=class_weight,
                                    random_state=random_state, solver=solver, max_iter=max_iter,
                                    multi_class=multi_class, verbose=verbose, warm_start=warm_start, n_jobs=n_jobs)
        # 使用同样的参数创建L2逻辑回归
        self.l2 = LogisticRegression(penalty='l2', dual=dual, tol=tol, C=C, fit_intercept=fit_intercept,
                                     intercept_scaling=intercept_scaling, class_weight=class_weight,
                                     random_state=random_state, solver=solver, max_iter=max_iter,
                                     multi_class=multi_class, verbose=verbose, warm_start=warm_start, n_jobs=n_jobs)

    def fit(self, X, y, sample_weight=None):
        # 训练L1逻辑回归
        super(LR, self).fit(X, y, sample_weight=sample_weight)
        self.coef_old_ = self.coef_.copy()
        # 训练L2逻辑回归
        self.l2.fit(X, y, sample_weight=sample_weight)

        cntOfRow, cntOfCol = self.coef_.shape
        # 权值系数矩阵的行数对应目标值的种类数目
        for i in range(cntOfRow):
            for j in range(cntOfCol):
                coef = self.coef_[i][j]
                # L1逻辑回归的权值系数不为0
                if coef != 0:
                    idx = [j]
                    # 对应在L2逻辑回归中的权值系数
                    coef1 = self.l2.coef_[i][j]
                    for k in range(cntOfCol):
                        coef2 = self.l2.coef_[i][k]
                        # 在L2逻辑回归中，权值系数之差小于设定的阈值，且在L1中对应的权值为0
                        if abs(coef1 - coef2) < self.threshold and j != k and self.coef_[i][k] == 0:
                            idx.append(k)
                    # 计算这一类特征的权值系数均值
                    mean = coef / len(idx)
                    self.coef_[i][idx] = mean
        return self


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
    y = [0] * len(group1) + [1] * len(group2)
    return X, y


def filter_variance(X, y, all_gene, args=None):
    variances = VarianceThreshold().fit(X, y)
    feature_gene = [gene for gene, variance in zip(all_gene, variances.variances_) if variance > args.threshold]
    return feature_gene


def Wrapper_LogisticRegression(X, y, all_gene, args=None):
    result = RFE(estimator=LogisticRegression(penalty=args.penalty, C=args.C), n_features_to_select=args.n_feature, step=args.step).fit(X, y)
    feature_gene = [gene for gene, variance in zip(all_gene, result.support_) if result == True]
    return feature_gene


def Embedded_L1_L2(X, y, all_gene, args=None):
    result = SelectFromModel(LR(threshold=0.5, C=0.1)).fit_transform(iris.data, iris.target)


def mean_decrease_impurity(X, y, all_gene, args=None):
    result = RandomForestRegressor(n_estimators=args.n_estimators, max_features=args.n_feature).fit(X, y)
    feature_gene = [gene for gene, importance in zip(all_gene, result.feature_importances_) if importance > args.threshold]
    return feature_gene


def Stability_selection(X, y, all_gene, args=None):
    result = RandomizedLasso(alpha=args.alpha).fit(X, y)
    feature_gene = [gene for gene, score in zip(all_gene, result.scores_) if score > args.threshold]
    return feature_gene


def lasso_regress(X, y, all_gene, args=None):
    result = Lasso(alpha=args.alpha).fit(X, y)
    print(result.coef_)
    feature_gene = [gene for gene, coef in zip(all_gene, result.coef_) if coef > args.threshold]
    return feature_gene


def logistic_regress(X, y, all_gene, args=None):
    result = LogisticRegression(penalty=args.penalty, C=args.C).fit(X, y)
    feature_gene = [gene for gene, coef in zip(all_gene, *result.coef_) if coef > args.threshold]
    return feature_gene


def feature_select(data, group1, group2, args=None):
    ## 选择特征基因，支持多种选择方法， 默认采用logistic回归拟合权重系数
    base_function = ["wilcox", "pearsonr", "fisher"]
    if args.feature_selection_method in base_function:
        num_group1 = len(group1); num_group2 = len(group2)
        pearsonr_target = [0] * num_group1 + [1] * num_group2
        result = {"gene": data.index, "a_vs_b": []}
        for x, y in zip(data[group1].as_matrix(), data[group2].as_matrix()):
            if args.feature_selection_method == "wilcox":
                z, p = ranksums(x, y)
            elif args.feature_selection_method == "pearsonr":
                cor, p = pearsonr(pearsonr_target, list(x) + list(y))
            elif args.feature_selection_method == "fisher":
                x_nonzero = np.count_nonzero(x); y_nonzero = np.count_nonzero(y)
                oddsratio, p = fisher_exact([[x_nonzero, num_group1 - x_nonzero], [y_nonzero, num_group2 - y_nonzero]])
            else:
                raise MethodException('must offer a right feature selection method')
            result["a_vs_b"].append(p)
        result = pd.DataFrame(result)
        result = result.loc[result["a_vs_b"] <= args.pval]["gene"]
        if len(result) == 0:
            raise MethodException("the number of feature gene is 0, please choose your fsm or other parameter")
        return list(result)
    else:
        sc = StandardScaler()
        genes = data.index
        X, y = get_train_test_data(data, group1, group2, genes)
        sc.fit(X)
        X_train_std = sc.transform(X)
        if args.feature_selection_method == "Lasso":
            feature_gene = lasso_regress(X_train_std, y, genes, args=args)
        elif args.feature_selection_method == "logistic":
            feature_gene = logistic_regress(X_train_std, y, genes, args=args)
        elif args.feature_selection_method == "RandomizedLasso":
            feature_gene = Stability_selection(X_train_std, y, genes, args=args)
        elif args.feature_selection_method == "RandomForest":
            feature_gene = mean_decrease_impurity(X_train_std, y, genes, args=args)
        elif args.feature_selection_method == "Wrapper":
            feature_gene = Wrapper_LogisticRegression(X_train_std, y, genes, args=args)
        elif args.feature_selection_method == "variance":
            feature_gene = filter_variance(X_train_std, y, genes, args=args)
        else:
            raise MethodException('must offer a right feature selection method')
        if len(feature_gene) == 0:
            raise MethodException("the number of feature gene is 0, please choose your fsm or other parameter")
        return feature_gene


def evaluate_model(data, group1, group2, feature_gene, args=None, name=None, method="LinearSVC", C=1, n_folds=5):
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

