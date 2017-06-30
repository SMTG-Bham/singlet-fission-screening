#!/usr/bin/env python

from __future__ import unicode_literals

import os
import sys
import tarfile
import shutil
import tempfile

from contextlib import contextmanager
from pymatgen.io.gaussian import GaussianInput, GaussianOutput

from tinydb import TinyDB


@contextmanager
def cd(run_path, cleanup=lambda: True):
    """
    Temporarily work in another directory, creating it if necessary.
    """
    home = os.getcwd()
    os.chdir(os.path.expanduser(run_path))
    try:
        yield
    finally:
        os.chdir(home)
        cleanup()


@contextmanager
def tempdir():
    """
    Temporarily work in temporary directory, deleting it aftewards.
    """
    dirpath = tempfile.mkdtemp()
    def cleanup():
        shutil.rmtree(dirpath)
    with cd(dirpath, cleanup):
        yield dirpath


def extract_data_from_tar_file(tar_file):
    with tarfile.open(tar_file, 'r:gz') as tar:
        tar.extractall()

    folder = tar_file.replace('.tar.gz', '')
    with cd(folder):
        tdout = GaussianOutput('td.log')
        td_exit = tdout.read_excitation_energies()
        td_triplet = [e for e in td_exit if 'triplet' in e[3].lower()][0][0]
        td_singlet = [e for e in td_exit if 'singlet' in e[3].lower()][0][0]

        tdaout = GaussianOutput('tda.log')
        tda_exit = tdaout.read_excitation_energies()
        tda_triplet = [e for e in tda_exit if 'triplet' in e[3].lower()][0][0]
        tda_singlet = [e for e in tda_exit if 'singlet' in e[3].lower()][0][0]

        nicssout = GaussianOutput('nics_singlet.log')
        # occasionally some jobs fail here
        if not nicssout.properly_terminated:
            return False
        nicss_mag = nicssout.read_magnetic_shielding()
        nicss_six_ring_above = (abs(nicss_mag[-8]['isotropic']) +
                                abs(nicss_mag[-6]['isotropic']))/2
        nicss_six_ring_below = (abs(nicss_mag[-7]['isotropic']) +
                                abs(nicss_mag[-5]['isotropic']))/2
        nicss_five_ring_above = (abs(nicss_mag[-4]['isotropic']) +
                                 abs(nicss_mag[-2]['isotropic']))/2
        nicss_five_ring_below = (abs(nicss_mag[-3]['isotropic']) +
                                 abs(nicss_mag[-1]['isotropic']))/2

        nicstout = GaussianOutput('nics_triplet.log')
        if not nicstout.properly_terminated:
            return False
        nicst_mag = nicstout.read_magnetic_shielding()
        nicst_six_ring_above = (abs(nicst_mag[-8]['isotropic']) +
                                abs(nicst_mag[-6]['isotropic']))/2
        nicst_six_ring_below = (abs(nicst_mag[-7]['isotropic']) +
                                abs(nicst_mag[-5]['isotropic']))/2
        nicst_five_ring_above = (abs(nicst_mag[-4]['isotropic']) +
                                 abs(nicst_mag[-2]['isotropic']))/2
        nicst_five_ring_below = (abs(nicst_mag[-3]['isotropic']) +
                                 abs(nicst_mag[-1]['isotropic']))/2

        data = {'td_singlet': td_singlet, 'td_triplet': td_triplet,
                'tda_singlet': tda_singlet, 'tda_triplet': tda_triplet,
                'nicss_six_ring_above': nicss_six_ring_above,
                'nicss_six_ring_below': nicss_six_ring_below,
                'nicss_five_ring_above': nicss_five_ring_above,
                'nicss_five_ring_below': nicss_five_ring_below,
                'nicst_six_ring_above': nicst_six_ring_above,
                'nicst_six_ring_below': nicst_six_ring_below,
                'nicst_five_ring_above': nicst_five_ring_above,
                'nicst_five_ring_below': nicst_five_ring_below}
    return data


data_to_write = []
db = TinyDB(os.path.join('..', 'data', 'structures.json'))
systems = list(db.all())
done = 0
for i, system in enumerate(systems):
    input_file = GaussianInput.from_dict(system['input'])
    directory = input_file.title
    tar_name = '{}.tar.gz'.format(directory)
    tar_file = os.path.abspath(os.path.join('..', 'data', 'calculations', tar_name))
    if os.path.isfile(tar_file):
        # extract the data in a temp directory to avoid clobbering any data
        with tempdir() as tmp_dir:
            shutil.copy(tar_file, tmp_dir)
            data = extract_data_from_tar_file(tar_name)
            if not data:
                print('{} did not finish correctly, skipping'.format(directory))
                continue
            data.update({'x_sub': system['x_sub'], 'y_sub': system['y_sub'],
                         'z_sub': system['z_sub'], 'nx': system['nx'],
                         'ny': system['ny'], 'title': system['title']})
            data_to_write.append(data)

    if i % 500 == 0:
        done += 5
        print('{}% completed'.format(done))

print('writing data')
db = TinyDB(os.path.join('..', 'data', 'calculated-data.json'))
db.insert_multiple(data_to_write)
