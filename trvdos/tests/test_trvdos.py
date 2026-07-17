"""
Unit and regression test for the trvdos package.
"""

# Import package, test suite, and other packages as needed
import sys
import trvdos as vd
import filecmp
import pickle
import MDAnalysis as mda
import pytest
from pathlib import Path
import numpy as np
import re
    
@pytest.fixture
def resprops_path(tmp_path):
    return tmp_path / "residueProperties.dat"
    
@pytest.fixture
def VACF_path(tmp_path):
    return tmp_path / "VACF.dat"

@pytest.fixture
def VDoS_path(tmp_path):
    return tmp_path / "VDoS.dat"


@pytest.fixture
def expected_singlestep():
    file_path = Path(__file__).parent / "expected_outputs" / "vdos_1step.pkl"
    with open(file_path, 'rb') as file:
        return pickle.load(file)

@pytest.fixture
def expected_fulltraj():
    file_path = Path(__file__).parent / "expected_outputs" / "vdos_fulltraj.pkl"
    with open(file_path, 'rb') as file:
        return pickle.load(file)


@pytest.fixture
def singlestep():
    TOPOL = Path(__file__).parent / "inputs/run-NVE.tpr"
    TRAJ =  Path(__file__).parent / "inputs/run-NVE_trunc.trr"
    u = mda.Universe(TOPOL,TRAJ)
    sel = u.select_atoms("resname POPC")
    vdos = vd.vdos(sel,10)
    vdos.single_frame(0, u.trajectory.time)
    return vdos

@pytest.fixture
def fulltraj(resprops_path, VACF_path, VDoS_path):
    TOPOL = Path(__file__).parent / "inputs/run-NVE.tpr"
    TRAJ =  Path(__file__).parent / "inputs/run-NVE_trunc.trr"
    u = mda.Universe(TOPOL,TRAJ)
    sel = u.select_atoms("resname POPC")
    vdos = vd.vdos(sel,10)
    tStep = 0
    for ts in (u.trajectory):
        vdos.single_frame(tStep,ts.time)
        tStep += 1
    return vdos

@pytest.fixture
def fulltraj_fileoutput(fulltraj, resprops_path, VACF_path, VDoS_path):
    # output files to temp path for later comparison
    fulltraj.copyResidueList()
    fulltraj.postProcess(fulltraj.residueListCopy, mode = "all")
    fulltraj.outputGeometry(resprops_path,fulltraj.residueListCopy)
    fulltraj.outputVACF(VACF_path)
    fulltraj.outputVDoS(VDoS_path)
    return fulltraj

@pytest.fixture
def expected_resprops_path():
    return Path(__file__).parent / "expected_outputs/residueProperties.dat"

@pytest.fixture
def expected_VACF_path():
    return Path(__file__).parent / "expected_outputs/VACF.dat"

@pytest.fixture
def expected_VDoS_path():
    return Path(__file__).parent / "expected_outputs/VDoS.dat"
    
    
def test_trvdos_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "trvdos" in sys.modules

def test_1step_totVACF(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('totVACF'), singlestep.totVACF)

def test_1step_totVDoS(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('totVDOS'), singlestep.totVDoS)
    
    
def test_1step_trVACF(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('trVACF'), singlestep.trVACF)

def test_1step_trVDoS(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('trVDOS'), singlestep.trVDoS)
    
    
def test_1step_rotVACF(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('rotVACF'), singlestep.rotVACF)

def test_1step_rotVDoS(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('rotVDOS'), singlestep.rotVDoS)
    
    
def test_1step_rotBondVACF(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('rotVACF'), singlestep.rotVACF)

def test_1step_rotBondVDoS(expected_singlestep, singlestep):
    assert np.array_equal(expected_singlestep.get('rotBondVDOS'), singlestep.rotBondVDoS)
    
    
def test_fulltraj_totVACF(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('totVACF'), fulltraj.totVACF)

def test_fulltraj_totVDoS(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('totVDOS'), fulltraj.totVDoS)
    
    
def test_fulltraj_trVACF(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('trVACF'), fulltraj.trVACF)

def test_fulltraj_trVDoS(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('trVDOS'), fulltraj.trVDoS)
    
    
def test_fulltraj_rotVACF(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('rotVACF'), fulltraj.rotVACF)

def test_fulltraj_rotVDoS(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('rotVDOS'), fulltraj.rotVDoS)
    
    
def test_fulltraj_rotBondVACF(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('rotVACF'), fulltraj.rotVACF)

def test_fulltraj_rotBondVDoS(expected_fulltraj, fulltraj):
    assert np.array_equal(expected_fulltraj.get('rotBondVDOS'), fulltraj.rotBondVDoS)

    
def test_residue_properties_output(expected_resprops_path, resprops_path, fulltraj_fileoutput):
    import re

    with open(expected_resprops_path) as f:
        expected = np.array([
            float(x)
            for x in re.findall(
                r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
                f.read(),
            )
        ])

    with open(resprops_path) as f:
        actual = np.array([
            float(x)
            for x in re.findall(
                r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?",
                f.read(),
            )
        ])

    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-10,
        atol=1e-12,
    )
    
def test_VACF_output(expected_VACF_path, VACF_path, fulltraj_fileoutput):
    expected = np.loadtxt(expected_VACF_path, comments="#")
    actual = np.loadtxt(VACF_path, comments="#")

    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-10,
        atol=1e-12,
    )
    
def test_VDoS_output(VDoS_path, expected_VDoS_path, fulltraj_fileoutput):
    expected = np.loadtxt(expected_VDoS_path, comments="#")
    actual = np.loadtxt(VDoS_path, comments="#")

    np.testing.assert_allclose(
        actual,
        expected,
        rtol=1e-6,
        atol=1e-3
    )