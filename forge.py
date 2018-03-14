#!/usr/bin/python3


import h5py
import click
import os


@click.group()
def main():
    pass


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-n', '--n-sims', type=int)
def ssplit(filename, n_sims):
    """Single run: Split a runfile created by plan.py into n-sims approx. equally large runfiles"""
    print(filename, n_sims)
    print(os.path.splitext(filename))
    name, extension = os.path.splitext(filename)
    f = h5py.File(filename, 'r')
    fs = [h5py.File(name + str(x+1) + extension, 'w') for x in range(n_sims)]

    for pf in fs:
        pf.create_group('Protocol')
        f.copy('Parameters', pf)
        pass

    i = 0
    for p in f['Protocol']:
        f['Protocol'].copy(p, fs[i]['Protocol'])
        i += 1
        i %= n_sims


@main.command()
@click.argument('filenames', nargs=-1)
@click.option('-o', '--outfile', type=click.Path())
def sjoin(filenames, outfile):
    """Single run: Join a number of outputfiles into one. Assumes they used the same systemfile"""
    print(filenames)
    f = h5py.File(outfile, 'w')
    f.create_group('Protocol')
    f.create_group('Record')
    fs = [h5py.File(x, 'r') for x in filenames]

    fs[0].copy('Parameters', f)
    fs[0].copy('Drugs', f)
    fs[0].copy('Gene', f)
    fs[0].copy('Rates', f)

    # joins lists of all protocolnames into one
    # we need this such that the records can be renamed correctly
    # as the name of a record is it's index in the protocol list
    protocolnames = sorted([z for y in [x['Protocol'] for x in fs] for z in y])
    print(protocolnames)
    print(len(protocolnames))

    for i, n in enumerate(protocolnames):
        ppath = 'Protocol/' + n
        for pf in fs:
            if ppath in pf:
                rpath = 'Record/0/' + str(list(pf['Protocol']).index(n))
                pf.copy(rpath, f['Record'], str(i))
                pf.copy(ppath, f['Protocol'], n)



if __name__ == "__main__":
    main()
