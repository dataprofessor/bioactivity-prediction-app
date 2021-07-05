import subprocess
import datetime
import os
import OrganaziedModleGenerator as om
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from chembl_webresource_client.new_client import new_client


protein = input("Please Enter a protin\n")  # acetylcholinesterase

def createFolder(name):

	now = datetime.datetime.now()
	time = now.strftime("%Y-%m-%d_%H:%M:%S")
	folderName = name + "-" + time
	os.mkdir(folderName)

	return folderName

def part1(protein, pathToSave):
    # Target search for a virus
    target = new_client.target
    target_query = target.search(protein)
    targets = pd.DataFrame.from_dict(target_query)
    selected_target = targets.target_chembl_id[0]
    activity = new_client.activity
    res = activity.filter(target_chembl_id=selected_target).filter(standard_type="IC50")
    df = pd.DataFrame.from_dict(res)
    df.to_csv(pathToSave+'/01_bioactivity_data_raw.csv', index=False)
    df2 = df[df.standard_value.notna()]
    df2 = df2[df.canonical_smiles.notna()]
    df2_nr = df2.drop_duplicates(['canonical_smiles'])
    selection = ['molecule_chembl_id', 'canonical_smiles', 'standard_value']
    df3 = df2_nr[selection]
    df3.to_csv(pathToSave+'/02_bioactivity_data_preprocessed.csv', index=False)
    df4 = pd.read_csv(pathToSave+'/02_bioactivity_data_preprocessed.csv')

    # classification
    bioactivity_threshold = []

    for i in df4.standard_value:
        if float(i) >= 10000:
            bioactivity_threshold.append("inactive")
        elif float(i) <= 1000:
            bioactivity_threshold.append("active")
        else:
            bioactivity_threshold.append("intermediate")

    bioactivity_class = pd.Series(bioactivity_threshold, name='class')
    df5 = pd.concat([df4, bioactivity_class], axis=1)
    df5.to_csv(pathToSave+'/03_bioactivity_data_curated.csv', index=False)

# Inspired by: https://codeocean.com/explore/capsules?query=tag:data-curation

def lipinski(smiles, verbose=False):

    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:

        desc_MolWt = Descriptors.MolWt(mol)
        desc_MolLogP = Descriptors.MolLogP(mol)
        desc_NumHDonors = Lipinski.NumHDonors(mol)
        desc_NumHAcceptors = Lipinski.NumHAcceptors(mol)

        row = np.array([desc_MolWt,
                        desc_MolLogP,
                        desc_NumHDonors,
                        desc_NumHAcceptors])

        if i == 0:
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i + 1

    columnNames = ["MW", "LogP", "NumHDonors", "NumHAcceptors"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors

def pIC50(input):
    pIC50 = []

    for i in input['standard_value_norm']:
        molar = i*(10**-9) # Converts nM to M
        pIC50.append(-np.log10(molar))

    input['pIC50'] = pIC50
    x = input.drop('standard_value_norm', 1)

    return x

def norm_value(input):
    norm = []

    for i in input['standard_value']:
        if i > 100000000:
          i = 100000000
        norm.append(i)

    input['standard_value_norm'] = norm
    x = input.drop('standard_value', 1)

    return x


def part2(pathToSave):
    df = pd.read_csv(pathToSave+'/03_bioactivity_data_curated.csv')
    df_no_smiles = df.drop(columns='canonical_smiles')
    smiles = []

    for i in df.canonical_smiles.tolist():
        cpd = str(i).split('.')
        cpd_longest = max(cpd, key=len)
        smiles.append(cpd_longest)

    smiles = pd.Series(smiles, name='canonical_smiles')

    df_clean_smiles = pd.concat([df_no_smiles, smiles], axis=1)

    df_lipinski = lipinski(df_clean_smiles.canonical_smiles)
    df_combined = pd.concat([df, df_lipinski], axis=1)

    df_combined.standard_value.describe()

    df_norm = norm_value(df_combined)
    df_norm.standard_value_norm.describe()
    df_final = pIC50(df_norm)
    df_final.pIC50.describe()
    df_final.to_csv(pathToSave+'/04_bioactivity_data_3class_pIC50.csv')
    df_2class = df_final[df_final['class'] != 'intermediate']
    df_2class.to_csv(pathToSave+'/05_bioactivity_data_2class_pIC50.csv')

def part3(pathToSave):
    df3 = pd.read_csv(pathToSave+'/04_bioactivity_data_3class_pIC50.csv')
    selection = ['canonical_smiles', 'molecule_chembl_id']
    df3_selection = df3[selection]
    df3_selection.to_csv('molecule.smi', sep='\t', index=False, header=False)
    #bash padel.sh
    print("In Bash")
    bashCommand = "bash padel.sh"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    # change the new file's location
    bashCommand = "mv descriptors_output.csv " + pathToSave + "/descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    df3_X = pd.read_csv(pathToSave+'/descriptors_output.csv')
    df3_X = df3_X.drop(columns=['Name'])
    df3_Y = df3['pIC50']
    dataset3 = pd.concat([df3_X,df3_Y], axis=1)
    dataset3.to_csv(pathToSave+'/06_bioactivity_data_3class_pIC50_pubchem_fp.csv', index=False)

def copyEssentialFiles(path):

    copyCommand = "cp " + path + "/descriptors_output.csv bioactivity-prediction-app-main"
    process = subprocess.Popen(copyCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    copyCommand = "cp " + path + "/descriptor_list.csv bioactivity-prediction-app-main"
    process = subprocess.Popen(copyCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    copyCommand = "cp " + path + "/modle.pkl bioactivity-prediction-app-main"
    process = subprocess.Popen(copyCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()


pathToSave = createFolder(protein)
part1(protein, pathToSave)
print("Part 2")
part2(pathToSave)
print("Part 3")
part3(pathToSave)
print("STEP 2 - Generate the modle")
om.run(pathToSave)
copyEssentialFiles(pathToSave)
