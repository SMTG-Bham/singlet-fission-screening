#!/usr/bin/env python

import os
import copy
import random
import itertools

from pymatgen import Composition
from pymatgen.io.gaussian import GaussianInput
from pymatgen.core.structure import Molecule


gin = GaussianInput.from_file('input-template.com')
gin.route_parameters['integral'] = '(acc2e=12)'  # hack to get round pmg bug

# define the substituents
subs = [('nitro', 'nitro'),
        ('amine', 'amine'),
        ('cyano', 'cyano'),
        ('hydroxyl', 'hydroxyl'),
        ('fluoro', 'fluoro'),
        ('chloro', Molecule(['X', 'Cl'], [[0., 0., 0.], [0., 0., 1.11]])),
        # ('bromo', Composition('Br')),
        ('cf3', Molecule(['X', 'C', 'F', 'F', 'F'],
                         [[-2.09, 0., -0.20], [-0.32, 0., -0.20],
                          [0.26, 0.83, -1.64], [0.26, 0.83, 1.24],
                          [0.26, -1.66, -0.20]]))]

# define the x, y, z substituent positions as {name: site_ids}, where name is
# consistent with Pakapol naming convention and site_id is a list of zero based
# array indexes
sub_sites = {'x': [48, 49], 'y': [50, 51], 'z': [52, 53]}
sub_sites_thiol = {'x': [47, 49], 'y': [30, 31], 'z': [32, 33]}

# these lists give the nitrogen positions as (name, site_ids), where name is
# consistent with Pakapol naming convention, and site_id is a list of zero based
# array indexes
nx_sites = [(1, [21, 28]),
            (2, [22, 27]),
            (3, [23, 26]),
            (4, [24, 25]),
            [None]]
ny_sites = [(1, [19, 30]),
            (2, [18, 32]),
            (3, [17, 33]),
            (4, [16, 34]),
            (5, [15, 35]),
            [None]]
nx_sites_thiol = [(1, [16, 23]),
                  (2, [17, 22]),
                  (3, [18, 21]),
                  (4, [19, 20]),
                  [None]]

# define the rings used to calculate NICS
six_mem_a = [0, 1, 2, 3, 4, 5]
six_mem_b = [1, 2, 6, 7, 9, 9]
five_mem_a = [0, 1, 9, 13, 12]
five_mem_b = [2, 3, 10, 11, 6]


def save_input_file(ginp, nx=None, ny=None, nz=None, x_sub=None, y_sub=None,
                    z_sub=None, prefix='ciba', folder='structure-files'):
    nx_str = nx if nx else ""
    ny_str = ny if ny else ""
    x_str = x_sub if x_sub else ""
    y_str = y_sub if y_sub else ""
    z_str = z_sub if z_sub else ""
    ginp.title = '{}_nx-{}_ny-{}_x-{}_y-{}_z-{}'.format(prefix, nx_str, ny_str,
                                                        x_str, y_str, z_str)
    ginp.write_file(os.path.join(folder, '{}.com'.format(ginp.title)),
                    cart_coords=True)


##############################################
# generate unsubstituted pyridine structures #
##############################################

for nx, ny in itertools.product(nx_sites, ny_sites):
    mol = copy.deepcopy(gin)

    # don't bother with unsubtituted structure
    if not nx[0] or not ny[0]:
        continue

    if nx[0]:
        mol.molecule[nx[1][0]]._species = Composition('N')
        mol.molecule[nx[1][1]]._species = Composition('N')
    if ny[0]:
        mol.molecule[ny[1][0]]._species = Composition('N')
        mol.molecule[ny[1][1]]._species = Composition('N')

    save_input_file(mol, nx=nx[0], ny=ny[0])

############################################
# generate substituted pyridine structures #
############################################

# randomly choose from nitrogen sites and substitutional sites
# one thing to consider is that nx positions 2 and 3 are incompatible
# with substitutional sites z and y, respectively and ny
# position 3 is incompatible with substitional site x
#
# the other important note is that we need to substitute the atoms from the
# highest site_id downwards, otherwise the order and index of the atoms will
# change
generated = 0
while generated < 100:
    mol = copy.deepcopy(gin)

    # a some Nones to the lists to water the # substituents down a little
    x_sub = random.choice(subs + [[None]]*4)
    y_sub = random.choice(subs + [[None]]*4)
    z_sub = random.choice(subs + [[None]]*4)
    nx = random.choice(nx_sites + [[None]]*3)
    ny = random.choice(ny_sites + [[None]]*3)

    # if there are any clashes or no substituents then skip the structure
    if ((x_sub[0] and ny[0] == 3) or (y_sub[0] and nx[0] == 2) or
            (z_sub[0] and nx[0] == 3)):
        continue
    elif not x_sub[0] and not y_sub[0] and not z_sub[0]:
        continue

    # these substituions happen in place therefore don't care about order
    if nx[0]:
        mol.molecule[nx[1][0]]._species = Composition('N')
        mol.molecule[nx[1][1]]._species = Composition('N')
    if ny[0]:
        mol.molecule[ny[1][0]]._species = Composition('N')
        mol.molecule[ny[1][1]]._species = Composition('N')

    # make a list of the substitutions and sort from highest to lowest site_id
    sub_list = []
    if x_sub[0]:
        sub_list += zip(sub_sites['x'], [x_sub[1]]*2)
    if y_sub[0]:
        sub_list += zip(sub_sites['y'], [y_sub[1]]*2)
    if z_sub[0]:
        sub_list += zip(sub_sites['z'], [z_sub[1]]*2)
    sub_list = sorted(sub_list, reverse=True, key=lambda x: x[0])

    # now we can make substitutions in correct order
    for site, sub in sub_list:
        # do this to get round really weird pmg bug, where a bromo group can't
        # be subtituted into the structure more than twice, whereas a fluoro and
        # chloro are fine (issue #687 in pymatgen github)
        # actually, had to scrap the work around, it doesn't work.
        # bromine is dead...
        mol.molecule.substitute(site, sub)

    save_input_file(mol, nx=nx[0], ny=ny[0], x_sub=x_sub[0], y_sub=y_sub[0],
                    z_sub=z_sub[0])
    generated += 1

#############################################
# generate substituted thiophene structures #
#############################################

# same as before, including substituent clashing and site_id ordering caveats
# but don't have a ny nitrogen substituent this time

gin = GaussianInput.from_file('input-template-thiol.com')
gin.route_parameters['integral'] = '(acc2e=12)'  # hack to get round pmg bug

generated = 0
while generated < 20:
    mol = copy.deepcopy(gin)

    # a some Nones to the lists to water the # substituents down a little
    x_sub = random.choice(subs + [[None]]*4)
    y_sub = random.choice(subs + [[None]]*4)
    z_sub = random.choice(subs + [[None]]*4)
    nx = random.choice(nx_sites_thiol + [[None]]*3)

    # if there are any clashes or no substituents then skip the structure
    if ((y_sub[0] and nx[0] == 2) or (z_sub[0] and nx[0] == 3)):
        continue
    elif not x_sub[0] and not y_sub[0] and not z_sub[0]:
        continue

    if nx[0]:
        mol.molecule[nx[1][0]]._species = Composition('N')
        mol.molecule[nx[1][1]]._species = Composition('N')

    # make a list of the substitutions and sort from highest to lowest site_id
    sub_list = []
    if x_sub[0]:
        sub_list += zip(sub_sites_thiol['x'], [x_sub[1]]*2)
    if y_sub[0]:
        sub_list += zip(sub_sites_thiol['y'], [y_sub[1]]*2)
    if z_sub[0]:
        sub_list += zip(sub_sites_thiol['z'], [z_sub[1]]*2)
    sub_list = sorted(sub_list, reverse=True, key=lambda x: x[0])

    for site, sub in sub_list:
        mol.molecule.substitute(site, sub)

    save_input_file(mol, nx=nx[0], x_sub=x_sub[0], y_sub=y_sub[0],
                    z_sub=z_sub[0], prefix='ciba_thiol')
    generated += 1
