# qsm_runner.py
import numpy as np
import PyTLidar.treeqsm as qsm
from PyTLidar import Utils

def run_qsm(file, ipd=0.05, pdmin=0.03, pdmax=0.12, brad1=0.06, brad2=0.13, plot=0):
    """
    Run QSM on a LAS file and return tree metrics.
    """
    # Load and normalize point cloud
    P = Utils.load_point_cloud(file)
    P = P - np.mean(P, axis=0)

    # Define inputs
    Inputs = qsm.define_input(P, 1, 1, 1)[0]
    Inputs['PatchDiam1'] = [ipd]
    Inputs['PatchDiam2Min'] = [pdmin]
    Inputs['PatchDiam2Max'] = [pdmax]
    Inputs['BallRad1'] = [brad1]
    Inputs['BallRad2'] = [brad2]
    Inputs['plot'] = plot

    # QSM pipeline
    cover1 = qsm.cover_sets(P, Inputs)
    cover1, Base, Forb = qsm.tree_sets(P, cover1, Inputs)
    segment1 = qsm.segments(cover1, Base, Forb)
    segment1 = qsm.correct_segments(P, cover1, segment1, Inputs, 0, 1, 1)
    RS = qsm.relative_size(P, cover1, segment1)

    cover2 = qsm.cover_sets(P, Inputs, RS)
    cover2, Base, Forb = qsm.tree_sets(P, cover2, Inputs, segment1)
    segment2 = qsm.segments(cover2, Base, Forb)
    segment2 = qsm.correct_segments(P, cover2, segment2, Inputs, 1, 1, 0)

    # Cylinders + branches + treedata
    cylinder = qsm.cylinders(P, cover2, segment2, Inputs)
    branch = qsm.branches(cylinder)
    treedata, _ = qsm.tree_data(cylinder, branch, P, Inputs, iter=1)

    return treedata.to_dict(orient="list")  # return as dict (easy for R)
