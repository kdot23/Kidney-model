__author__ = 'zhuoshuli'

import functions
from DistributionGenerator import DistributionGenerator
import util
import copy

import random
import numpy as np
import pdb

class SaidmanCompatibleGenerator:
    def __init__(self):

        self.Pr_FEMALE = 0.4090
        self.Pr_SPOUSAL_DONOR = 0.4897
        self.Pr_LOW_PRA = 0.7019
        self.Pr_MED_PRA = 0.2

        self.Pr_LOW_PRA_INCOMPATIBILITY = 0.05
        self.Pr_MED_PRA_INCOMPATIBILITY = 0.45
        self.Pr_HIGH_PRA_INCOMPATIBILITY = 0.90

        self.Pr_SPOUSAL_PRA_COMPATIBILITY = 0.75

        self.Pr_PATIENT_TYPE_O = 0.4814
        self.Pr_PATIENT_TYPE_A = 0.3373
        self.Pr_PATIENT_TYPE_B = 0.1428

        self.Pr_DONOR_TYPE_O = 0.4814
        self.Pr_DONOR_TYPE_A = 0.3373
        self.Pr_DONOR_TYPE_B = 0.1428
        self.Pr_DONOR_TYPE_AB = 1 - 0.4814 - 0.3373 - 0.1428

        self.LKDPI_mu = 37.1506024096
        self.LKDPI_std = 22.2170610307

        self.LKDPI_mu_off = 40.10
        self.LKDPI_std_off = 20.49

        self.LKDPI_mu_1 = 35.23
        self.LKDPI_std_1 = 19.45

        self.LKDPI_mu_2 = 37.76
        self.LKDPI_std_2 = 19.63

        self.LKDPI_mu_3 = 41.89
        self.LKDPI_std_3 = 20.33

        self.LKDPI_mu_4 = 53.24
        self.LKDPI_std_4 = 20.02

        self.type_O_given_AA = 0.51
        self.type_A_given_AA = 0.257
        self.type_B_given_AA = 0.19
        self.type_AB_given_AA = 0.043

        self.donorAA = 592 / 5538.0

        self.AA_given_type_O = self.type_O_given_AA * self.donorAA / self.Pr_DONOR_TYPE_O
        self.AA_given_type_A = self.type_A_given_AA * self.donorAA / self.Pr_DONOR_TYPE_A
        self.AA_given_type_B = self.type_B_given_AA * self.donorAA / self.Pr_DONOR_TYPE_B
        self.AA_given_type_AB = self.type_AB_given_AA * self.donorAA / self.Pr_DONOR_TYPE_AB

    def drawPatientBloodType(self):

        # 'O': 0, 'A': 1, 'B': 2, 'AB': 3

        r = random.random()
        if (r <= self.Pr_PATIENT_TYPE_O): return 0
        if (r <= self.Pr_PATIENT_TYPE_O + self.Pr_PATIENT_TYPE_A): return 1
        if (r <= self.Pr_PATIENT_TYPE_O + self.Pr_PATIENT_TYPE_A + self.Pr_PATIENT_TYPE_B): return 2
        return 3

    def drawDonorBloodType(self):
        # 'O': 0, 'A': 1, 'B': 2, 'AB': 3
        r = random.random()
        if (r <= self.Pr_DONOR_TYPE_O): return 0
        if (r <= self.Pr_DONOR_TYPE_O + self.Pr_DONOR_TYPE_A): return 1
        if (r <= self.Pr_DONOR_TYPE_O + self.Pr_DONOR_TYPE_A + self.Pr_DONOR_TYPE_B): return 2
        return 3

    def drawAfricanAmerican(self, donor_blood_type):

        # 'O': 0, 'A': 1, 'B': 2, 'AB': 3
        if donor_blood_type == 0 :
            x = np.random.choice([0, 1], p=[1 - self.AA_given_type_O, self.AA_given_type_O])
        if donor_blood_type == 1:
            x = np.random.choice([0, 1], p=[1 - self.AA_given_type_A, self.AA_given_type_A])
        if donor_blood_type == 2:
            x = np.random.choice([0, 1], p=[1 - self.AA_given_type_B, self.AA_given_type_B])
        if donor_blood_type == 3:
            x = np.random.choice([0, 1], p=[1 - self.AA_given_type_AB, self.AA_given_type_AB])

        return x.astype(int)

    def isPatientFemale(self):
        return random.random() <= self.Pr_FEMALE

    def isDonorSpouse(self):
        return random.random() <= self.Pr_SPOUSAL_DONOR

    def isPositiveCrossmatch(self, pr_PraIncompatibility):
        return random.random() <= pr_PraIncompatibility

    def generatePraIncompatibility(self):
        r = random.random()
        if (r <= self.Pr_LOW_PRA):
            pr_PraIncompatiblity = self.Pr_LOW_PRA_INCOMPATIBILITY
        elif (r <= self.Pr_LOW_PRA + self.Pr_MED_PRA):
            pr_PraIncompatiblity = self.Pr_MED_PRA_INCOMPATIBILITY
        else:
            pr_PraIncompatiblity = self.Pr_HIGH_PRA_INCOMPATIBILITY

        return pr_PraIncompatiblity

    def generateCompPraIncompatibility(self, isWifePatient, pr_PraIncompatiblity):
        if (~isWifePatient):
            return pr_PraIncompatiblity
        else:
            return 1.0 - self.Pr_SPOUSAL_PRA_COMPATIBILITY * (1.0 - pr_PraIncompatiblity)

    def generateCompatibility(self, bloodTypeDonor, bloodTypePatient, patientCPRA):
        return functions.are_blood_compatible(bloodTypeDonor, bloodTypePatient) \
               and (not self.isPositiveCrossmatch(patientCPRA))

    def generateLKDPI(self):
        # x = 100
        # while x > 93 or x < -8:
        x = random.gauss(self.LKDPI_mu, self.LKDPI_std)
        return int(x)

    def generateLKDPIdonor(self, m, lkdpi):
        if lkdpi <= 20:
            x = np.random.normal(self.LKDPI_mu_1, self.LKDPI_std_1, m)
        elif lkdpi > 20 and lkdpi <= 40:
            x = np.random.normal(self.LKDPI_mu_2, self.LKDPI_std_2, m)
        elif lkdpi > 40 and lkdpi <= 60:
            x = np.random.normal(self.LKDPI_mu_3, self.LKDPI_std_3, m)
        else:
            x = np.random.normal(self.LKDPI_mu_4, self.LKDPI_std_4, m)

        return x.astype(int)

    # def generateLKDPIOff(self):
    #     # x = 200
    #     # while x > 105 or x < -20:
    #     x = random.gauss(self.LKDPI_mu_off, self.LKDPI_std_off)
    #     return int(x)

    # def generateLKDPIMat(self, m, n):
    #     x = np.random.normal(self.LKDPI_mu_off, self.LKDPI_std_off, (m, n))
    #     return x.astype(int)


class SaidmanPair():
    def __init__(self):
        self.generator = SaidmanCompatibleGenerator()
        self.generatePair()

    def generatePair(self):

        self.bloodTypePatient = self.generator.drawPatientBloodType()
        self.bloodTypeDonor = self.generator.drawDonorBloodType()
        self.pr_PraIncompatiblity = self.generator.generatePraIncompatibility()

        if (self.generator.isPatientFemale()):
            self.patientFemale = True
        else:
            self.patientFemale = False

        if (self.generator.isDonorSpouse()):
            self.donorSpouse = True
        else:
            self.donorSpouse = False

        self.isWifePatient = self.patientFemale and self.donorSpouse
        self.patientCPRA = self.generator.generateCompPraIncompatibility(self.isWifePatient, self.pr_PraIncompatiblity)

        self.compatible = functions.are_blood_compatible(self.bloodTypeDonor, self.bloodTypePatient) \
                          and (not self.generator.isPositiveCrossmatch(self.patientCPRA))

        self.LKDPI = self.generator.generateLKDPI()

class SaidmanLKDPIPool:

    def __init__(self, size):
        self.generatePool(size)

    def generatePool(self, size):

        self.allPairs = []
        self.compatiblePairs = []
        self.incompatiblePairs = []

        self.allID = []
        self.compatibleID = []
        self.incompatibleID = []

        for i in range(size):
            pair = None
            pair = SaidmanPair()

            self.allPairs.append(pair)
            self.allID.append(i)
            if pair.compatible:
                self.compatiblePairs.append(pair)
                self.compatibleID.append(i)
            else:
                self.incompatiblePairs.append(pair)
                self.incompatibleID.append(i)

class BJCPair():
    def __init__(self, life_time=0, pair_index=0):
        self.life_time = life_time
        self.pair_index = pair_index
        self.generator = DistributionGenerator()
        self.saidman = SaidmanCompatibleGenerator()
        self.generatePair()
        self.get_lkdpi()

    def generatePair(self):
        self.bloodTypePatient = self.saidman.drawPatientBloodType()
        self.bloodTypeDonor = self.saidman.drawDonorBloodType()
        self.donor_afam = self.saidman.drawAfricanAmerican(self.bloodTypeDonor)

        self.donor_age = self.generator.gen_donor_age()
        self.donor_egfr = self.generator.gen_donor_egfr(self.donor_age)
        self.donor_sbp = self.generator.gen_donor_sbp(self.donor_egfr)
        self.donor_sex = self.generator.gen_donor_sex()
        # self.donor_afam = self.generator.gen_donor_afam()

        self.donor_cig_use = self.generator.gen_donor_cig_use()
        # 1 is male
        self.rec_sex = self.generator.gen_rec_sex()
        self.donor_weight = self.generator.gen_donor_weight(self.donor_sex)
        self.rec_weight = self.generator.gen_rec_weight(self.rec_sex)
        self.donor_bmi = self.generator.gen_donor_bmi(self.donor_weight)

        self.pr_PraIncompatiblity = self.saidman.generatePraIncompatibility()
        
        if (self.rec_sex == 0):
            self.patientFemale = True
        else:
            self.patientFemale = False
        
        
        if (self.saidman.isDonorSpouse()):
            self.donorSpouse = True
        else:
            self.donorSpouse = False

        self.isWifePatient = self.patientFemale and self.donorSpouse
        self.patientCPRA = self.saidman.generateCompPraIncompatibility(self.isWifePatient, self.pr_PraIncompatiblity)

        self.compatible = functions.are_blood_compatible(self.bloodTypeDonor, self.bloodTypePatient) \
                          and (not self.saidman.isPositiveCrossmatch(self.patientCPRA))

        self.donor_rec_related = self.generator.gen_donor_rec_related()
        self.donor_rec_abo_comp = functions.are_blood_compatible(self.bloodTypeDonor, self.bloodTypePatient) #self.generator.gen_not_donor_rec_abo_comp()
        # self.donor_rec_HLA_B_mis = self.generator.gen_donor_rec_HLA_B_mis(self.donor_rec_related)
        # self.donor_rec_HLA_DR_mis = self.generator.gen_donor_rec_HLA_DR_mis(self.donor_rec_related)
        mis_B_DR = self.generator.gen_donor_rec_HLA_mis(self.donor_rec_related)
        self.donor_rec_HLA_B_mis = mis_B_DR[0]
        self.donor_rec_HLA_DR_mis = mis_B_DR[1]

        self.donor_rec_weight_ratio = util.calculate_weight_ratio(self.donor_weight, self.rec_weight)


    def get_lkdpi(self):
        self.LKDPI = util.calculate_lkdpi(self.donor_age, self.donor_afam, self.donor_bmi, self.donor_cig_use,
                                          self.donor_sex, self.rec_sex,
                                          self.donor_sbp, self.donor_rec_abo_comp,
                                          self.donor_rec_related, self.donor_egfr,
                                          self.donor_rec_HLA_B_mis, self.donor_rec_HLA_DR_mis,
                                          self.donor_rec_weight_ratio)

        self.survival = util.calculate_survival(self.LKDPI)

class BJCpool():
    def __init__(self, size):
        self.generatePool(size)

    def generatePool(self, size):

        self.allPairs = []
        self.compatiblePairs = []
        self.incompatiblePairs = []

        self.allID = []
        self.compatibleID = []
        self.incompatibleID = []

        self.num_incomp_O_type = 0
        self.num_comp_O_type = 0

        for i in range(size):
            pair = None
            pair = BJCPair()

            self.allPairs.append(copy.deepcopy(pair))
            self.allID.append(i)
            if pair.compatible:
                self.compatiblePairs.append(copy.deepcopy(pair))
                self.compatibleID.append(i)
                if pair.bloodTypePatient == 0:
                    self.num_comp_O_type += 1
            else:
                self.incompatiblePairs.append(copy.deepcopy(pair))
                self.incompatibleID.append(i)
                if pair.bloodTypePatient == 0:
                    self.num_incomp_O_type += 1


class BJCSensitivityPool:

    def __init__(self, compSize, incompSize):
        self.generatePool(compSize, incompSize)

    def generatePool(self, compSize, incompSize):
        self.allPairs = []
        self.compatiblePairs = []
        self.incompatiblePairs = []

        self.compatibleIDs = []
        self.incompatibleIDs = []

        self.num_incomp_O_type = 0
        self.num_comp_O_type = 0

        index = 0

        count = 0
        while count < compSize:
            pair = None
            pair = BJCPair(pair_index=index)

            if pair.compatible:
                self.allPairs.append(copy.deepcopy(pair))
                self.compatiblePairs.append(copy.deepcopy(pair))
                self.compatibleIDs.append(index)
                if pair.bloodTypePatient == 0:
                    self.num_comp_O_type += 1
                index += 1
                count += 1

        count = 0
        while count < incompSize:
            pair = None
            pair = BJCPair(pair_index=index)

            if not pair.compatible:
                self.allPairs.append(copy.deepcopy(pair))
                self.incompatiblePairs.append(copy.deepcopy(pair))
                self.incompatibleIDs.append(index)
                if pair.bloodTypePatient == 0:
                    self.num_incomp_O_type += 1
                index += 1
                count += 1


class BJCSensitivityPoolRandomize:

    def __init__(self, compSize, incompSize):
        self.generatePool(compSize, incompSize)

    def generatePool(self, compSize, incompSize):

        self.allPairs = [None] * (compSize + incompSize)
        self.compatiblePairs = []
        self.incompatiblePairs = []

        tmp_ID = list(np.random.permutation(compSize+incompSize))
        self.allID = range(compSize+incompSize)
        self.compatibleID = []
        self.incompatibleID = []

        index = 0
        count = 0
        while count < compSize:
            pair = BJCPair()
            if pair.compatible:
                real_id = tmp_ID[index]
                self.allPairs[real_id] = copy.deepcopy(pair)
                self.compatiblePairs.append(copy.deepcopy(pair))
                self.compatibleID.append(real_id)
                index += 1
                count += 1

        count = 0
        while count < incompSize:
            pair = BJCPair()
            if not pair.compatible:
                real_id = tmp_ID[index]
                self.allPairs[real_id] = copy.deepcopy(pair)
                self.incompatiblePairs.append(copy.deepcopy(pair))
                self.incompatibleID.append(real_id)
                index += 1
                count += 1

