import multiprocessing
from pathlib import Path

import numpy as np
import prody as pd
from Bio.PDB.ResidueDepth import get_surface, residue_depth, ca_depth, min_dist

import helper


def get_residue_depths(pdb_file):
    """
    Get residue depths

    Parameters
    ----------
    pdb_file

    Returns
    -------
    dict of depth _ ca/cb/mean
    """
    structure, residues, _, _, _ = helper.read_pdb(pdb_file)
    surface = get_surface(structure)
    data = {"depth_mean": np.array([residue_depth(residue, surface) for residue in residues]),
            "depth_cb": np.array([min_dist(helper.get_beta_coordinates(residue), surface) for residue in residues]),
            "depth_ca": np.array([ca_depth(residue, surface) for residue in residues])}
    return data


def get_fluctuations(protein: pd.AtomGroup, n_modes: int = 50):
    """
    Get atom fluctuations using anisotropic and Gaussian network models with n_modes modes.

    Parameters
    ----------
    protein
    n_modes

    Returns
    -------
    dict of anm_ca, anm_cb, gnm_ca, gnm_cb
    """
    data = {}
    beta_indices = helper.get_beta_indices(protein)
    alpha_indices = helper.get_alpha_indices(protein)
    data["anm_cb"] = get_anm_fluctuations(protein[beta_indices], n_modes)
    data["gnm_cb"] = get_gnm_fluctuations(protein[beta_indices], n_modes)
    data["anm_ca"] = get_anm_fluctuations(protein[alpha_indices], n_modes)
    data["gnm_ca"] = get_gnm_fluctuations(protein[alpha_indices], n_modes)
    return data


def get_anm_fluctuations(protein: pd.AtomGroup, n_modes: int = 50):
    """
    Get atom fluctuations using an Anisotropic network model with n_modes modes.
    """
    protein_anm, _ = pd.calcANM(protein, n_modes=n_modes, selstr='all')
    return pd.calcSqFlucts(protein_anm)


def get_gnm_fluctuations(protein: pd.AtomGroup, n_modes: int = 50):
    """
    Get atom fluctuations using a Gaussian network model with n_modes modes.
    """
    protein_gnm, _ = pd.calcGNM(protein, n_modes=n_modes, selstr='all')
    return pd.calcSqFlucts(protein_gnm)


def get_features_multiple(pdb_files, dssp_dir, num_threads=20, only_dssp=True, force_overwrite=False):
    """
    Extract features for a list of pdb_files in parallel

    Parameters
    ----------
    pdb_files
    dssp_dir
        directory to store tmp dssp files
    num_threads
    only_dssp
        extract only dssp features (use if not interested in features)
    force_overwrite
        force rerun DSSP

    Returns
    -------
    List of feature dicts (same order as pdb_files)
    """
    num_threads = min(len(pdb_files), num_threads)
    with multiprocessing.Pool(processes=num_threads) as pool:
        return pool.starmap(get_features, [(pdb_file, dssp_dir, only_dssp, force_overwrite) for pdb_file in pdb_files])


def get_features(pdb_file: str, dssp_dir: str, only_dssp=True, force_overwrite=True):
    """
    Extract features from a pdb_file

    Parameters
    ----------
    pdb_file
    dssp_dir
        directory to store tmp dssp files
    only_dssp
        extract only dssp features (use if not interested in features)
    force_overwrite
        force rerun DSSP

    Returns
    -------
    dict of features
    """
    pdb_file = str(pdb_file)
    _, name, _ = helper.get_file_parts(pdb_file)
    protein = pd.parsePDB(pdb_file).select("protein")
    # if Path(pdb_file).suffix != ".pdb":
    pdb_file = str(Path(dssp_dir) / f"{name}.pdb")
    pd.writePDB(pdb_file, protein)
    protein = pd.parsePDB(pdb_file)
    dssp_file = Path(dssp_dir) / f"{name}.dssp"
    if force_overwrite or not dssp_file.exists():
        dssp_file = pd.execDSSP(str(pdb_file), outputname=name, outputdir=str(dssp_dir))
    protein = pd.parseDSSP(dssp=str(dssp_file), ag=protein, parseall=True)
    data = get_dssp_features(protein)
    if only_dssp:
        return data
    else:
        data = {**data, **get_fluctuations(protein)}
        data = {**data, **get_residue_depths(pdb_file)}
        # data = {**data, **get_electrostatics(protein, pdb_file, es_dir=dssp_dir)}
        return data


def get_dssp_features(protein_dssp):
    """
    Extracts DSSP features (assumes DSSP is run already)

    Parameters
    ----------
    protein_dssp
        protein on which execDSSP has been called
    Returns
    -------
    dict of secondary,
    dssp_
        NH_O_1_index, NH_O_1_energy
            hydrogen bonds; e.g. -3,-1.4 means: if this residue is residue i then N-H of I is h-bonded to C=O of I-3 with an
            electrostatic H-bond energy of -1.4 kcal/mol. There are two columns for each type of H-bond, to allow for bifurcated H-bonds.
        NH_O_2_index, NH_O_2_energy
        O_NH_1_index, O_NH_1_energy
        O_NH_2_index, O_NH_2_energy
        acc
            number of water molecules in contact with this residue *10. or residue water exposed surface in Angstrom^2.
        alpha
            virtual torsion angle (dihedral angle) defined by the four C?? atoms of residues I-1,I,I+1,I+2. Used to define chirality.
        kappa
            virtual bond angle (bend angle) defined by the three C?? atoms of residues I-2,I,I+2. Used to define bend (structure code ???S???).
        phi
            IUPAC peptide backbone torsion angles.
        psi
            IUPAC peptide backbone torsion angles.
        tco
            cosine of angle between C=O of residue I and C=O of residue I-1. For ??-helices, TCO is near +1, for ??-sheets TCO is near -1.
    Ignores:
    dssp_bp1, dssp_bp2, and dssp_sheet_label: residue number of first and second bridge partner followed by one letter sheet label
    """
    dssp_ignore = ["dssp_bp1", "dssp_bp2", "dssp_sheet_label", "dssp_resnum"]
    dssp_labels = [label for label in protein_dssp.getDataLabels() if label.startswith("dssp") and label not in dssp_ignore]
    data = {}
    alpha_indices = helper.get_alpha_indices(protein_dssp)
    indices = [protein_dssp[x].getData("dssp_resnum") for x in alpha_indices]
    for label in dssp_labels:
        label_to_index = {i - 1: protein_dssp[x].getData(label) for i, x in zip(indices, alpha_indices)}
        data[f"{label}"] = np.array([label_to_index[i] if i in label_to_index else 0 for i in range(len(alpha_indices))])
    data["secondary"] = protein_dssp.getData("secondary")[alpha_indices]
    return data


