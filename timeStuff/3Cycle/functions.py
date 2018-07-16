""" Functions used in prepareMatching.py """

import numpy as np
import pdb
import pandas as pd


def calculate_LKDPI(x):
	""" Calculates a LKDPI score when given vector x. Note that the parameters must be in the order below.
	 * 0. Donor age
	 * 1. Donor sex 									
	 * 2. Recipeint sex									
	 * 3. Donor eGFR
	 * 4. Donor SBP
	 * 5. Donor BMI
	 * 6. Donor is African-American						
	 * 7. Donor history of cigarette use				
	 * 8. Donor and recipient are *RELATED*			
	 * 9. Donor and recipient are ABO *COMPATIBLE*	
	 * 10. Donor and recipient weight ratio
	 * 11. Donor and recipient HLA-B mismatches
	 * 12. Donor and recipient HLA-DR mismatches
	 * 13. Donor ID
	 * 14. Recipient ID """

	age = x[0]
	sbp = x[4]
	bmi = x[5]
    
	#Fix age
	age = age - 50 if age > 50 else 0

	#Fix bmi
	if(x[6] == 1):									#if African-American
		bmi += 22.34
	if(x[7] == 1):									#if history of cigarette use
		bmi += 14.33
		    
	#Fix sbp
	if(bool(x[1] == 1) & bool(x[2] == 1)):			#if both male
		sbp -= 21.68
	if(x[9] == 0):									#if ABO incompatible
		sbp += 27.30
	if(x[8] == 0):			  						#if unrelated
		sbp -= 10.61

  	#Calculate LKDPI
	# lkdpi = -11.30 + 1.85 * age - 0.381 * x[4] + 1.17 * bmi + 0.44 * sbp + 8.57 * x[11] + 8.26 * x[12] - 50.87 * x[10]
	# find a bug:
	lkdpi = -11.30 + 1.85 * age - 0.381 * x[3] + 1.17 * bmi + 0.44 * sbp + 8.57 * x[11] + 8.26 * x[12] - 50.87 * x[10]
	return int(np.ceil(lkdpi))

def calculate_survival(y):
	""" Calculates the expected survival for a vector, y, of LKDPI scores. """

	return 14.78*np.exp(-0.01239*y['LKDPI'])

def return_donor_col():
	""" Returns columns of the original data frame that are donor specific"""

	return [3, 4, 5, 6, 7, 8, 9, 11, 13, 15, 17, 32]

def return_recipient_col():
	""" Returns columns of the original data frame that are donor specific"""
	return [0, 1, 2, 12, 14, 16]

def are_blood_compatible(donor_blood, recipient_blood):
	""" Returns 1 if donor and recipient are ABO compatible, 0 otherwise. Note that the input for this function
		is NOT a vector. """
        
	if(donor_blood == 0):
		return 1
	if(donor_blood == recipient_blood):
		return 1
	if(recipient_blood == 3):
		return 1
	return 0

def calculate_mimatches(donor, recipient):
	recipient_list = [set(string.split("/")) for string in recipient]
	donor_list =  [set(string.split("/")) for string in donor]
	result = []
	for i in range(len(recipient_list)):
		recipient = recipient_list[i]
		if '' in recipient:
			recipient.remove('')
		donor = donor_list[i]
		if '' in donor:
			donor.remove('')
		difference = len(donor.difference(recipient))
		bool_one = difference <= 1 and len(donor) == 1
		bool_two = len(recipient) == 1 and difference == 0
		if  bool_one and (not bool_two):
			difference = difference + 1
		result.append(difference)

	return result

def fill_dataset(x):
	copy_frame = x.copy()
	""" Returns dataset with Donor/Rec weight ratio, ABO Compatibility, HLA-B Mimatches, and HLA-DR Mismatches calculated. """

	#Determine ABO compatibility
	copy_frame['Donor/Rec ABO *COMPATIBLE*'] = map(lambda y,z: are_blood_compatible(y,z), copy_frame['Donor Blood Type'].values, copy_frame['Recipient Blood Type'].values)
	
	#Calculate weight ratio and return minimum
	copy_frame['Donor/Rec Weight Ratio'] = np.minimum(copy_frame['Donor Weight'] / copy_frame['Recipient Weight'], 0.9)

	#Calculate HLA mismatches
	copy_frame['Donor/Rec HLA-B Mismatches'] = calculate_mimatches(copy_frame['Donor HLA-B'], copy_frame['Recipient HLA-B'])
	copy_frame['Donor/Rec HLA-DR Mismatches'] = calculate_mimatches(copy_frame['Donor HLA-DR'], copy_frame['Recipient HLA-DR'])

	# copy_frame.drop(['Donor Blood Type', 'Recipient Blood Type', 'Donor Weight', 'Recipient Weight', 'Donor HLA-B', 'Recipient HLA-B', 'Donor HLA-DR', 'Recipient HLA-DR' ], axis = 1, inplace = True)
	# copy_frame.drop(['Donor Blood Type', 'Recipient Blood Type', 'Donor HLA-B', 'Recipient HLA-B', 'Donor HLA-DR', 'Recipient HLA-DR' ], axis = 1, inplace = True)
	copy_frame.drop(['Donor HLA-B', 'Recipient HLA-B', 'Donor HLA-DR',
					 'Recipient HLA-DR'], axis=1, inplace=True)

	return copy_frame
	

# donors = ['a/b', 'b/a','a/','a/b','a/b','a/','a/','a/b','a/b', 'a/b','a/']
# recipients = ['a/b', 'b/a','a/','a/c','c/b','a/b','b/a','a/','b/', 'c/d','b/']

# sd = pd.Series(donors)
# sr = pd.Series(recipients)
# result = calculate_mimatches(sd, sr)
# print result

def return_file_paths(condition):
	""" Returns list of file paths for 'Original', 'Optimal', and 'Counterfactual' files based on condition.

		condition: 'None', 'Adjusted', 'Blood', 'Both'
	"""
	if condition == 'None':
		return ["~/Dropbox/working/Sofie/Results/Original Pairs/OriginalPairsCalculated.csv", "~/Dropbox/working/Sofie/Results/No Conditions/FinalResultsNoConditions.csv", "~/Dropbox/working/Sofie/Results/Counterfactual Pairs/CounterfactualPairsCalculated.csv"]
	if condition == 'Adjusted':
		return ["~/Dropbox/working/Sofie/Results/Original Pairs/OriginalPairsCalculated.csv", "~/Dropbox/working/Sofie/Results/Must Improve LKDPI/FinalResultsAdjusted.csv", "~/Dropbox/working/Sofie/Results/Counterfactual Pairs/CounterfactualPairsCalculated.csv"]
	if condition == 'Blood':
		return ["~/Dropbox/working/Sofie/Results/Original Pairs/OriginalPairsCalculated.csv", "~/Dropbox/working/Sofie/Results/No ABO Incompatible/FinalResultsABO.csv", "~/Dropbox/working/Sofie/Results/Counterfactual Pairs/CounterfactualPairsCalculated.csv"]
	if condition == 'Both':
		return ["~/Dropbox/working/Sofie/Results/Original Pairs/OriginalPairsCalculated.csv", "~/Dropbox/working/Sofie/Results/Both Conditions/FinalResultsBothConditions.csv", "~/Dropbox/working/Sofie/Results/Counterfactual Pairs/CounterfactualPairsCalculated.csv"]
	return ["**You did not specify a condition**"]