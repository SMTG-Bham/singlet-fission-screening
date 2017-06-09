#!/usr/bin/env python


import os
import copy
import itertools

from pymatgen import Composition
from pymatgen.io.gaussian import GaussianInput
from pymatgen.core.structure import Molecule

from tinydb import TinyDB


# define the substituents
subs = [('nitro', 'nitro'),
        ('amine', 'amine'),
        ('cyano', 'cyano'),
        ('hydroxyl', 'hydroxyl'),
        ('fluoro', 'fluoro'),
        ('chloro', Molecule(['X', 'Cl'], [[0., 0., 0.], [0., 0., 1.11]])),
        ('bromo', Molecule(['X', 'Br'], [[0., 0., 0.], [0., 0., 1.11]])),
        [None]]

# define the x, y, z substituent positions as {name: site_ids}, where name is
# consistent with Pakapol naming convention and site_id is a list of zero based
# array indexes
sub_sites = {'x': [48, 49], 'y': [50, 51], 'z': [52, 53]}
sub_sites_thiol = {'x': [47, 49], 'y': [30, 31], 'z': [32, 33]}

# these lists give the nitrogen positions as (name, N site_ids, H site_ids),
# where name is consistent with Pakapol naming convention, and site_id is a
# list of zero based array indexes
nx_sites = [(1, [21, 28], [43, 44]),
            (2, [22, 27], [52, 53]),
            (3, [23, 26], [50, 51]),
            (4, [24, 25], [45, 42]),
            [None]]
ny_sites = [(1, [19, 30], [36, 40]),
            (2, [18, 32], [41, 37]),
            (3, [17, 33], [48, 49]),
            (4, [16, 34], [46, 47]),
            (5, [15, 35], [38, 39]),
            [None]]
nx_sites_thiol = [(1, [16, 23], [27, 28]),
                  (2, [17, 22], [32, 33]),
                  (3, [18, 21], [30, 31]),
                  (4, [19, 20], [26, 29]),
                  [None]]

# define the rings used to calculate NICS
six_mem_a = [0, 1, 2, 3, 4, 5]
six_mem_b = [1, 2, 6, 7, 9, 9]
five_mem_a = [0, 1, 9, 13, 12]
five_mem_b = [2, 3, 10, 11, 6]

data_to_write = []
index = 1


def cache_input_file(ginp, index, nx=None, ny=None, x_sub=None, y_sub=None,
                     z_sub=None, prefix='ciba'):
    nx_str = nx if nx else ""
    ny_str = ny if ny else ""
    x_str = x_sub if x_sub else ""
    y_str = y_sub if y_sub else ""
    z_str = z_sub if z_sub else ""
    ginp.title = '{}_nx-{}_ny-{}_x-{}_y-{}_z-{}'.format(prefix, nx_str, ny_str,
                                                        x_str, y_str, z_str)
    data_to_write.append({'input': ginp.as_dict(), 'nx': nx, 'ny': ny,
                          'x_sub': x_sub, 'y_sub': y_sub, 'z_sub': z_sub,
                          'title': ginp.title, 'index': index})


############################################
# generate substituted pyridine structures #
############################################

# need to consider that nx positions 2 and 3 are incompatible
# with substitutional sites z and y, respectively and ny
# position 3 is incompatible with substitional site x
#
# the other important note is that we need to substitute the atoms from the
# highest site_id downwards, otherwise the order and index of the atoms will
# change

gin = GaussianInput.from_file(os.path.join('..', 'templates',
                                           'input-template.com'))
gin.route_parameters['integral'] = '(acc2e=12)'  # hack to get round pmg bug
    yout
print "generating pyridine substituted structures..."
for nx, ny, x_sub, y_sub, z_sub in itertools.product(nx_sites, ny_sites,
                                                     subs, subs, subs):
    mol = copy.deepcopy(gin)

    # if there are any clashes or no substituents then skip the structure
    if ((x_sub[0] and ny[0] == 3) or (y_sub[0] and nx[0] == 2) or
            (z_sub[0] and nx[0] == 3)):
        continue
    elif not x_sub[0] and not y_sub[0] and not z_sub[0]:
        continue

    # make a list of the substitutions and sort from highest to lowest site_id
    # this list includes the H atoms next to the N we introduce in to the rings
    sub_list = []

    # these substituions happen in place therefore don't care about order
    if nx[0]:
        mol.molecule[nx[1][0]]._species = Composition('N')
        mol.molecule[nx[1][1]]._species = Composition('N')
        sub_list += zip(nx[2], [None]*2)
    if ny[0]:
        mol.molecule[ny[1][0]]._species = Composition('N')
        mol.molecule[ny[1][1]]._species = Composition('N')
        sub_list += zip(ny[2], [None]*2)

    if x_sub[0]:
        sub_list += zip(sub_sites['x'], [x_sub[1]]*2)
    if y_sub[0]:
        sub_list += zip(sub_sites['y'], [y_sub[1]]*2)
    if z_sub[0]:
        sub_list += zip(sub_sites['z'], [z_sub[1]]*2)
    sub_list = sorted(sub_list, reverse=True, key=lambda x: x[0])

    # now we can make substitutions in correct order
    for site, sub in sub_list:
        if sub:
            mol.molecule.substitute(site, sub)
        else:
            # hyrdogen atom to remove
            gin.molecule.remove_sites([site])

    cache_input_file(mol, index, nx=nx[0], ny=ny[0], x_sub=x_sub[0],
                     y_sub=y_sub[0], z_sub=z_sub[0])
    index += 1

#############################################
# generate substituted thiophene structures #
#############################################

# same as before, including substituent clashing and site_id ordering caveats
# but don't have a ny nitrogen substituent this time

gin = GaussianInput.from_file(os.path.join('..', 'templates',
                                           'input-template-thiol.com'))
gin.route_parameters['integral'] = '(acc2e=12)'  # hack to get round pmg bug

print "generating thiophene substituted structures..."
for nx, x_sub, y_sub, z_sub in itertools.product(nx_sites, subs, subs, subs):
    mol = copy.deepcopy(gin)

    # if there are any clashes or no substituents then skip the structure
    if ((y_sub[0] and nx[0] == 2) or (z_sub[0] and nx[0] == 3)):
        continue
    elif not x_sub[0] and not y_sub[0] and not z_sub[0]:
        continue

    # make a list of the substitutions and sort from highest to lowest site_id
    sub_list = []

    if nx[0]:
        mol.molecule[nx[1][0]]._species = Composition('N')
        mol.molecule[nx[1][1]]._species = Composition('N')
        sub_list += zip(nx[2], [None]*2)

    if x_sub[0]:
        sub_list += zip(sub_sites_thiol['x'], [x_sub[1]]*2)
    if y_sub[0]:
        sub_list += zip(sub_sites_thiol['y'], [y_sub[1]]*2)
    if z_sub[0]:
        sub_list += zip(sub_sites_thiol['z'], [z_sub[1]]*2)
    sub_list = sorted(sub_list, reverse=True, key=lambda x: x[0])

    for site, sub in sub_list:
        if sub:
            mol.molecule.substitute(site, sub)
        else:
            # hyrdogen atom to remove
            gin.molecule.remove_sites([site])

    cache_input_file(mol, index, nx=nx[0], x_sub=x_sub[0], y_sub=y_sub[0],
                     z_sub=z_sub[0], prefix='ciba_thiol')
    index += 1


db = TinyDB(os.path.join('..', 'data', 'structures.json'))
db.insert_multiple(data_to_write)
