from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, AllChem

def analyze_compound(compound_name):
    """
    Analyze compound properties using RDKit.
    """
    try:
        sdf_filename = f"{compound_name}_3D.sdf"
        mol = Chem.SDMolSupplier(sdf_filename)[0]  # Load the first molecule from the SDF file

        if not mol:
            raise ValueError("Molecule could not be parsed from SDF.")

        # Generate 3D conformer
        AllChem.EmbedMultipleConfs(mol, numConfs=1)

        # Calculate properties
        mol_weight = Descriptors.MolWt(mol)
        smiles = Chem.MolToSmiles(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        hydrogen_bond_donors = Descriptors.NumHDonors(mol)
        hydrogen_bond_acceptors = Descriptors.NumHAcceptors(mol)
        polar_surface_area = Descriptors.TPSA(mol)

        # Generate 2D structure
        img = Draw.MolToImage(mol)
        img.save(f"{compound_name}_2D.png")
        print(f"2D structure saved as {compound_name}_2D.png")

        return {
            "Molecular Weight": mol_weight,
            "SMILES": smiles,
            "Rotatable Bonds": rotatable_bonds,
            "Hydrogen Bond Donors": hydrogen_bond_donors,
            "Hydrogen Bond Acceptors": hydrogen_bond_acceptors,
            "Polar Surface Area": polar_surface_area,
        }
    except Exception as e:
        print(f"Error analyzing compound: {e}")
    return None
