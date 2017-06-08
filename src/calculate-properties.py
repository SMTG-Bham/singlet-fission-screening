#!/usr/bin/env python

import os
import sys
import tarfile

import numpy as np
import scipy.linalg
import logging

from contextlib import contextmanager

from pymatgen.core.periodic_table import DummySpecie
from pymatgen.io.gaussian import GaussianInput, GaussianOutput

from tinydb import TinyDB, Query


@contextmanager
def cd(run_path):
    """
    Temporarily work in another directory, creating it if necessary.
    """
    home = os.getcwd()
    if not os.path.isdir(run_path):
        os.makedirs(run_path)
    os.chdir(run_path)
    yield
    os.chdir(home)


def add_dummy_atoms(mol, ring_atom_ids):
    data = np.array([mol[i].coords for i in ring_atom_ids])
    centre = np.average(data, axis=0)

    # fit the atomic positions to a plane using least squares method, taken
    # from: https://gist.github.com/amroamroamro/1db8d69b4b65e8bc66a6#file-curve_fitting-py-L27
    A = np.c_[data[:, 0], data[:, 1], np.ones(data.shape[0])]
    C, _, _, _ = scipy.linalg.lstsq(A, data[:, 2])

    # calculate the normal unit vector
    normal = np.array(np.cross([1, 0, C[0]], [0, 1, C[1]]))
    unit = normal/scipy.linalg.norm(normal)

    # add the dummy atoms above and below the ring
    mol.append(DummySpecie('X-Bq'), centre + unit)
    mol.append(DummySpecie('X-Bq'), centre - unit)
    return mol


# define calculation settings
functional = 'b3lyp'
dieze_tag = '#p'
basis_set = '6-311++G(d,p)'
link0 = {'%oldchk': 'chkpt.chk', '%mem': '58GB', '%nprocshared': '24'}
td_params = {'integral': '(acc2e=12)', 'td': '(50-50)', 'guess': 'read'}
tda_params = {'integral': '(acc2e=12)', 'tda': '(50-50)', 'guess': 'read'}
nmr_params = {'integral': '(acc2e=12)', 'guess': 'read', 'nmr': ''}

# define the rings used to calculate NICS
six_mem_a = [0, 1, 2, 3, 4, 5]
six_mem_b = [1, 2, 6, 7, 8, 9]
five_mem_a = [0, 1, 9, 13, 12]
five_mem_b = [2, 3, 10, 11, 6]

index = int(sys.argv[1])

# use absolute path so we don't loose track of the db when changing directory
db = TinyDB(os.path.abspath('structures.json'))
query = Query()
compound = db.get(query.index == index)
rin = GaussianInput.from_dict(compound['input'])
directory = rin.title


def calculate_properties():
    logging.info('started processing {}'.format(directory))

    with cd(directory):
        rin.write_file('relax.com', cart_coords=True)
        os.system('g09 < relax.com > relax.log')
        rout = GaussianOutput('relax.log')
        if not rout.properly_terminated:
            logging.error('{} relaxation did not terminate correctly'.format(
                          directory))
            return

        # do the TD-DFT calculation
        tdin = GaussianInput(rout.final_structure, charge=0, title=rin.title,
                             functional=functional, spin_multiplicity=1,
                             basis_set=basis_set, dieze_tag=dieze_tag,
                             link0_parameters=link0, route_parameters=td_params)
        tdin.write_file('td.com', cart_coords=True)
        os.system('g09 < td.com > td.log')
        tdout = GaussianOutput('td.log')
        if not tdout.properly_terminated:
            logging.error('{} TD-DFT did not terminate correctly'.format(
                          directory))
            return

        # do the TD-DFT calculation w. Tamm-Dancoff approx
        tdain = GaussianInput(rout.final_structure, charge=0, title=rin.title,
                              spin_multiplicity=1, functional=functional,
                              basis_set=basis_set, dieze_tag=dieze_tag,
                              link0_parameters=link0,
                              route_parameters=tda_params)
        tdain.write_file('tda.com', cart_coords=True)
        os.system('g09 < tda.com > td.log')
        tdaout = GaussianOutput('tda.log')
        if not tdaout.properly_terminated:
            logging.error('{} TDA-DFT did not terminate correctly'.format(
                          directory))
            return

        # add the dummy atoms for the NICS(1)_zz calculations
        mol_nics = rout.final_structure
        mol_nics = add_dummy_atoms(mol_nics, six_mem_a)
        mol_nics = add_dummy_atoms(mol_nics, six_mem_b)
        mol_nics = add_dummy_atoms(mol_nics, five_mem_a)
        mol_nics = add_dummy_atoms(mol_nics, five_mem_b)

        # run NICS on the ground state and triplet state
        nicssin = GaussianInput(mol_nics, charge=0, title=rin.title,
                                spin_multiplicity=1, functional=functional,
                                basis_set=basis_set, dieze_tag=dieze_tag,
                                link0_parameters=link0,
                                route_parameters=nmr_params)
        nicssin.write_file('nics_singlet.com', cart_coords=True)
        # work around as pymatgen does not allow Bq as an element
        os.system("sed -i 's/X-Bq0+/Bq/g' nics_singlet.com")
        os.system('g09 < nics_singlet.com > nics_singlet.log')
        nicssout = GaussianOutput('nics_singlet.log')
        # can't have NICS job completion check due to above mentioned pmg bug

        nicstin = GaussianInput(mol_nics, charge=0, title=rin.title,
                                spin_multiplicity=3, functional=functional,
                                basis_set=basis_set, dieze_tag=dieze_tag,
                                link0_parameters=link0,
                                route_parameters=nmr_params)
        nicstin.write_file('nics_triplet.com', cart_coords=True)
        os.system("sed -i 's/X-Bq0+/Bq/g' nics_triplet.com")
        os.system('g09 < nics_triplet.com > nics_triplet.log')
        nicstout = GaussianOutput('nics_triplet.log')
    logging.info('finished processing {}'.format(directory))

with cd('calculations'):
    logging.basicConfig(filename='calculations.log',level=logging.DEBUG)

    calculate_properties()
    os.remove(os.path.join(directory, 'chkpt.chk'))
    with tarfile.open("{}.tar.gz".format(directory), "w:gz") as tar:
        tar.add(directory)
    shutil.rmtree(directory)