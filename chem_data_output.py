# import rdkit
from rdkit import Chem
# Chemical conformer, descriptor module, draw module
from rdkit.Chem import Draw, Descriptors, AllChem

def analyze_compound(compound_name):
    try:
        sdf_filename = f"{compound_name}_3D.sdf"
        # Load molecule from SDF file
        mol = Chem.SDMolSupplier(sdf_filename)[0]

        if not mol:
            # Error handling
            raise ValueError("Molecule could not be parsed from SDF.")

        # 3D conformer, based from databases
        AllChem.EmbedMultipleConfs(mol, numConfs=1)

        # Obtain properties
        # Molecular weight
        mol_weight = Descriptors.MolWt(mol)
        # SMILES string
        smiles = Chem.MolToSmiles(mol)
        # Number of rotatable bonds
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        # Number of hydrogen bond
        hydrogen_bond_donors = Descriptors.NumHDonors(mol)
        # Number of hydrogen bond acceptors
        hydrogen_bond_acceptors = Descriptors.NumHAcceptors(mol)
        # Polar surface area
        polar_surface_area = Descriptors.TPSA(mol)

        # Generates 2D structure
        img = Draw.MolToImage(mol)
        # Image save name format
        img.save(f"{compound_name}_2D.png")
        # Provides feedback
        print(f"2D structure saved as {compound_name}_2D.png")

        # Outputs molecular information
        return {
            "Molecular Weight": mol_weight,
            "SMILES": smiles,
            "Rotatable Bonds": rotatable_bonds,
            "Hydrogen Bond Donors": hydrogen_bond_donors,
            "Hydrogen Bond Acceptors": hydrogen_bond_acceptors,
            "Polar Surface Area": polar_surface_area,
        }
    except Exception as e:
        # Error handling
        print(f"Error analyzing compound: {e}")
    #     Return nothing when error
    return None
