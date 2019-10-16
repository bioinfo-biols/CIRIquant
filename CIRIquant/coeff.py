#! /usr/bin/env python
# -*- encoding:utf-8 -*-
import logging
import numpy as np


LOGGER = logging.getLogger('CIRIquant')


def correction(sample_exp, sample_stat, rnaser_exp, rnaser_stat):
    """
    RNase R correction
    """
    LOGGER.info('RNase R treatment coefficient correction')

    header = []
    circ_exp = {}

    train_index = []
    overlap_index = []

    for i in sample_exp:
        if sample_exp[i]['bsj'] <= 1:
            continue
        if i not in rnaser_exp or rnaser_exp[i]['bsj'] <= 5:
            continue
        if float(sample_exp[i]['bsj']) / sample_stat[1] >= float(rnaser_exp[i]['bsj']) / rnaser_stat[1]:
            continue
        overlap_index.append(i)

        if sample_exp[i]['ratio'] == 1 or rnaser_exp[i]['ratio'] == 1:
            continue
        if sample_exp[i]['ratio'] >= rnaser_exp[i]['ratio']:
            continue
        train_index.append(i)
    if len(overlap_index) == 0 or len(train_index) == 0:
        LOGGER.warn('No enough overlap circRNAs for correction, skipping this step')
        return [], sample_exp

    # Overlap expression constant for BSJ
    total_median = np.median([sample_exp[i]['bsj'] for i in overlap_index])
    rnaser_median = np.median([rnaser_exp[i]['bsj'] for i in overlap_index])
    exp_constant = total_median / rnaser_median

    y_factor = [factor(sample_exp[i]) for i in train_index]
    y, y_mean = np.transpose(np.matrix(y_factor)), np.median(y_factor)

    x_factor = [factor(rnaser_exp[i]) for i in train_index]
    x, x_mean = np.transpose(np.matrix(x_factor)), np.median(x_factor)
    coef_mean = y_mean / x_mean

    LOGGER.info('Fitting Model')
    coef, intercept = fit_model(x, y)
    LOGGER.info('Generate prior distribution ..')
    gmm = prior_distribution(x, y)

    header += [
        'RNaseR_Reads: {}'.format(rnaser_stat[1]),
        'Amplification: {}'.format(exp_constant),
        'Coef: {}'.format(coef),
        'Intercept: {}'.format(intercept),
        'N: {}'.format(gmm.n_components),
        'W: {}'.format(','.join([str(i) for i in np.round(gmm.weights_, 4)])),
        'M: {}'.format(','.join([str(i[0]) for i in np.round(gmm.means_, 4)])),
        'SD: {}'.format(','.join([str(i[0][0]) for i in np.round(gmm.covariances_, 4)])),
    ]

    for i in rnaser_exp:
        if i in sample_exp and sample_exp[i]['bsj'] == 0 and rnaser_exp[i]['fsj'] != 0:
            corrected_bsj = coef_mean * sample_exp[i]['fsj'] * rnaser_exp[i]['bsj'] / rnaser_exp[i]['fsj']
            corrected_ratio = junc_ratio(corrected_bsj, sample_exp[i]['fsj'])
            circ_exp[i] = {
                'bsj': corrected_bsj, 'fsj': sample_exp[i]['fsj'], 'ratio': corrected_ratio,
                'rnaser_bsj': rnaser_exp[i]['bsj'], 'rnaser_fsj': rnaser_exp[i]['fsj']
            }
        elif i in sample_exp and sample_exp[i]['bsj'] != 0:
            circ_exp[i] = sample_exp[i]
            # corrected_bsj = coef_mean * sample_exp[i]['fsj'] * rnaser_exp[i]['bsj'] / rnaser_exp[i]['fsj']
            # corrected_ratio = junc_ratio(corrected_bsj, sample_exp[i]['fsj'])
            # circ_exp[i] = {
            #     'bsj': corrected_bsj, 'fsj': sample_exp[i]['fsj'], 'ratio': corrected_ratio,
            #     'rnaser_bsj': rnaser_exp[i]['bsj'], 'rnaser_fsj': rnaser_exp[i]['fsj']
            # }
        else:
            corrected_bsj = rnaser_exp[i]['bsj'] * exp_constant
            corrected_fsj = corrected_bsj / (rnaser_exp[i]['bsj'] / rnaser_exp[i]['fsj']) / coef_mean if rnaser_exp[i]['fsj'] != 0 else 0
            corrected_ratio = junc_ratio(corrected_bsj, corrected_fsj)
            circ_exp[i] = {
                'bsj': corrected_bsj, 'fsj': corrected_fsj, 'ratio': corrected_ratio,
                'rnaser_bsj': rnaser_exp[i]['bsj'], 'rnaser_fsj': rnaser_exp[i]['fsj'],
            }

    for i in sample_exp:
        if i in circ_exp:
            continue
        circ_exp[i] = sample_exp[i]

    return header, circ_exp


def junc_ratio(bsj, fsj):
    return 2.0 * bsj / (2.0 * bsj + fsj)


def factor(d):
    """
    Factor
    """
    return d['ratio'] / (1.0 - d['ratio'])


def fit_model(x, y):
    from sklearn import linear_model, model_selection
    import warnings
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        x_train, x_test, y_train, y_test = model_selection.train_test_split(x, y, test_size=0.4, random_state=0)
        clf = linear_model.LinearRegression(fit_intercept=True)
        clf.fit(x_train, y_train)

        coef = clf.coef_[0][0]
        try:
            intercept = clf.intercept_[0]
        except Exception:
            intercept = 0
        LOGGER.debug('Coefficient: {}'.format(coef))
        LOGGER.debug('Intercept: {}'.format(intercept))
        LOGGER.debug('R-square: {}'.format(clf.score(x_train, y_train)))
        LOGGER.debug('Var: {}'.format((np.array(y_test - clf.predict(x_test)) ** 2).sum()))
    return coef, intercept


def prior_distribution(x, y):
    from sklearn.mixture import GaussianMixture
    import warnings

    coeff = y / x
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=DeprecationWarning)
        models = [GaussianMixture(i).fit(coeff) for i in np.arange(1, 4)]
    AIC = [m.aic(coeff) for m in models]
    # BIC = [m.bic(coeff) for m in models]
    gmm = models[np.argmin(AIC)]
    return gmm
