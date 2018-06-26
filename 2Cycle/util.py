__author__ = 'zhuoshuli'
import numpy as np
import pdb

def checkMatch(antigen, antibody):
    if len(antibody) == 0:
        return True
    if len(antibody) == 0:
        return True

    for item in antibody:
        if (item in antigen):
            return False

    return True


def calculate_survival(y):
    """ Calculates the expected survival for a vector, y, of LKDPI scores. """

    return 14.78*np.exp(-0.01239*y)


def calculate_mimatches(donor, recipient):
    recipient = set(recipient.split("/"))
    if '' in recipient:
        recipient.remove('')
    donor = set(donor.split("/"))
    if '' in donor:
        donor.remove('')
    difference = len(donor.difference(recipient))
    bool_one = difference <= 1 and len(donor) == 1
    bool_two = len(recipient) == 1 and difference == 0
    if bool_one and (not bool_two):
        difference = difference + 1
    return difference

def calculate_lkdpi(age, African_American, bmi, history_of_cigarette_use, donor_sex, recepient_sex,
                    sbp, ABO_incompatible, bio_related, eGFR, HLA_B_mismatches, HLA_DR_mismatches, D_R_WR):
    # Fix age
    if age > 50:
        age -= 50
    else:
        age = 0

    # Fix bmi
    if African_American == 1:
        bmi += 22.34
    if history_of_cigarette_use == 1:
        bmi += 14.33

    # Fix sbp
    if donor_sex == 1 and recepient_sex == 1: # Both male
        sbp -= 21.68
    if ABO_incompatible == 0:
        sbp += 27.30
    if bio_related == 0:
        sbp -= 10.61

    # Calculate LKDPI
    lkdpi = -11.30 + 1.85 * age - 0.381 * eGFR + 1.17 * bmi + 0.44 * sbp + \
            8.57 * HLA_B_mismatches + 8.26 * HLA_DR_mismatches - 50.87 * np.min([D_R_WR, 0.9])
    if lkdpi < -77:
        lkdpi = -77
    if lkdpi > 110:
        lkdpi = 110

    return int(np.ceil(lkdpi))

def calculate_weight_ratio(donor_weight, rec_weight):

    return np.minimum(1.0 * donor_weight/rec_weight, 0.9)

def calculate_LKDPI(s):

    return - (np.log(s) - np.log(14.87)) * 1 / 0.01239

def safe_div(x,y):
    if y == 0:
        return 0
    return x / y

def get_utility(exchange, w_ti=None, w_ij=None):

    comp_utility = 0
    incomp_utility = 0

    if w_ti is not None:
        T = w_ti.shape[0]
        N = w_ti.shape[1]

        for t in range(T):
            for i in range(N):
                real_t = exchange.pool.compatiblePairs[t].pair_index
                if w_ti[t, i] == 1:
                    if i == 0:
                        comp_utility += exchange.pool.compatiblePairs[t].survival
                    else:
                        real_i = exchange.pool.compatiblePairs[i-1].pair_index
                        comp_utility += exchange.survival_matrix[real_t, real_i]
                        incomp_utility += exchange.survival_matrix[real_i, real_t]

    if w_ij is not None:
        N = w_ij.shape[0]

        for i in range(N):
            for j in range(N):
                if w_ij[i, j] == 1:
                    real_i = exchange.pool.compatiblePairs[i].pair_index
                    real_j = exchange.pool.compatiblePairs[j].pair_index
                    incomp_utility += exchange.survival_matrix[real_i, real_j]
                    incomp_utility += exchange.survival_matrix[real_j, real_i]

    return comp_utility, incomp_utility
