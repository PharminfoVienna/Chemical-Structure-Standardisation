#############################################################################################################
# This script is developed since 2017 at the university of Vienna in the Pharmacoinformatics research group
#
# Description: The script can be used within a KNIME Python node to standardise chemical structures
# 
#
# Authors:       Jennifer Hemmerich (jennifer.hemmerich[at]univie.ac.at)
#
# Copyright 2020 Jennifer Hemmerich and Pharmacoinformatics research Group
#
#
##############################################################################################################



try:
    from rdkit import Chem
    from rdkit.Chem.SaltRemover import SaltRemover
    from rdkit.Chem.MolStandardize import rdMolStandardize
except ImportError:
    print(
        "rdkit could not be imported, please make sure that you have python 3 with rdkit version 2019 or higher installed")

import pandas as pd


def is_nonorganic(fragment):
    """Return true if fragment contains at least one carbon atom.
    :param fragment: The fragment as an RDKit Mol object.
    """
    # adapted from MolVS functiopn is_organic!!
    # TODO: Consider a different definition?
    # Could allow only H, C, N, O, S, P, F, Cl, Br, I
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() == 6:
            return False
    return True


def contains_nonorg(fragment):
    # organic: H, C, N, O, P, S, F, Cl, Br, I
    for a in fragment.GetAtoms():
        if a.GetAtomicNum() not in [1, 6, 7, 8, 15, 16, 9, 17, 35, 53]:
            return "Yes"
    return "No"


r = SaltRemover()

molecule_column = input_table['Molecule']  # Input from KNIME table
stand_mol_list = []
errs = []
mixture = "No"

for index, input_cell in molecule_column.iteritems():  # iterate through molecule list
    mol = input_cell
    if mol is None:
        stand_mol_list.append(("Got empty molecule", index, mol, "No", None, None))
        continue
    try:
        mol = rdMolStandardize.MetalDisconnector().Disconnect(mol)  # Disconnect metals
    except ValueError as e:
        if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) > 1:
            mixture = "Yes"
        stand_mol_list.append(("Failed at disconnect", index, None, mixture, None, str(e)))
        continue

    mol = r.StripMol(mol)

    # Check if we have multiple fragments present

    if len(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)) > 1:
        mixture = "Yes"
    else:
        mixture = "No"

    # Standardize fragments separately

    for i, frag in enumerate(Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)):

        frag = r.StripMol(frag)
        if frag.GetNumAtoms() == 0:
            continue
        elif is_nonorganic(frag):
            continue
        else:
            nonorg = contains_nonorg(frag)

            try:
                frag = rdMolStandardize.Normalize(frag)
            except ValueError as e:
                stand_mol_list.append(("Failed at normalize", index, None, mixture, nonorg, str(e)))
                continue
            try:
                frag = rdMolStandardize.Uncharger().uncharge(frag)
            except ValueError as e:
                stand_mol_list.append(("Failed at neutralising", index, None, mixture, nonorg, str(e)))
                continue

            if flow_variables['stereo'] == "Remove":
                try:
                    Chem.RemoveStereochemistry(frag)
                except ValueError as e:
                    stand_mol_list.append(("Failed at stereochem remove", index, None, mixture, nonorg, str(e)))
                    continue

            if flow_variables['stereo'] == "Clean":
                ''' 
                From RDKit documentation:
                Does the CIP stereochemistry assignment
    
                for the molecule’s atoms (R/S) and double bond (Z/E). Chiral atoms will have a property ‘_CIPCode’ indicating their chiral code.
    
                ARGUMENTS:
    
                mol: the molecule to use
                cleanIt: (optional) if provided, atoms with a chiral specifier that aren’t actually chiral (e.g. atoms with duplicate substituents or only 2 substituents, etc.) will have their chiral code set to CHI_UNSPECIFIED. Bonds with STEREOCIS/STEREOTRANS specified that have duplicate substituents based upon the CIP atom ranks will be marked STEREONONE.
                force: (optional) causes the calculation to be repeated, even if it has already been done
                flagPossibleStereoCenters (optional) set the _ChiralityPossible property on atoms that are possible stereocenters
    
                '''
                try:
                    Chem.AssignStereochemistry(frag, force=True, cleanIt=True)
                except ValueError as e:
                    stand_mol_list.append(("Failed at stereochem", index, None, mixture, nonorg, str(e)))
                    continue

        stand_mol_list.append(("Standardized", index, frag, mixture, nonorg, None))

df1 = input_table.copy()  # fetch input data
df1['Molecule_index'] = df1.index.astype(str).values  # Save Index as str column for merge

df2 = pd.DataFrame(stand_mol_list)
df2.columns = ['Standardized_success', 'Molecule_index', 'Molecule_standardised', 'Mixture', 'Nonorganic',
               'RDkit_error']  # set column names
df2['Molecule_index'] = df2['Molecule_index'].astype(str).values  # convert unicode column to string for merge

all = pd.merge(df1, df2, on='Molecule_index', how='right', validate="one_to_many")  # merge dataframes
all.index = range(len(all['Molecule_index']))  # set new index to prevent duplicate row ids

output_table_1 = all
output_table_2 = pd.DataFrame()

if flow_variables['Keep_all'] == "No":
    output_table_1 = all[all['Standardized_success'] == "Standardized"]
    out1 = all[all['Standardized_success'] != "Standardized"]
else:
    out1 = pd.DataFrame()

if flow_variables['keep_mixtures'] == "No":
    output_table_1 = output_table_1[output_table_1['Mixture'] == "No"]
    out2 = all[all['Mixture'] != "No"]
else:
    out2 = pd.DataFrame()

if flow_variables["keep_nonorganic"] == "No":
    output_table_1 = output_table_1[output_table_1['Nonorganic'] == "No"]
    out3 = all[all['Nonorganic'] != "No"]
else:
    out3 = pd.DataFrame()

output_table_2 = pd.concat([out1, out2, out3])
output_table_2 = output_table_2.groupby(by=output_table_2.index).first()
