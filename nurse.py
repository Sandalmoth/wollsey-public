#!/usr/bin/python3


import h5py
import click
import re
import numpy as np
import sys


TEXT_ENCODING = 'utf-8'


def make_rectangular(ll, fill):
    maxlength = max([len(l) for l in ll])
    for l in ll:
        while(len(l) < maxlength):
            l.append(fill)


def encode_ll(ll):
    nll = []
    for l in ll:
        nll.append([bytes(x, TEXT_ENCODING) for x in l])
    return nll


def decode_ll(ll):
    nll = []
    for l in ll:
        nll.append([x.decode(TEXT_ENCODING) for x in l])
    return nll


def ndiffs(a, b):
    nd = 0
    for x, y in zip(a, b):
        if x != y:
            nd += 1
    return nd


@click.group()
def main():
    pass


@main.command()
@click.argument('filename', type=click.Path())
def init(filename):
    """Run first to create a hdf5 file."""
    f = h5py.File(filename, 'w')
    f.create_group('Drugs')
    f.create_group('Gene')
    f.create_group('Rates')
    f.create_group('Cells')


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-i', '--min-cells', type=int)
@click.option('-a', '--max-cells', type=int)
@click.option('-z', '--population_size', type=int)
@click.option('-k', '--spring_force', type=float)
@click.option('-m', '--max-iters', type=int)
@click.option('-p', '--mutation_probability', type=float)
@click.option('-s', '--stop_condition', type=click.Choice(['MAX_ITERS', 'WT_HALF', 'WT_NULL', 'DEFAULT']), default='DEFAULT')
def setp(filename, min_cells, max_cells, max_iters, population_size, spring_force, mutation_probability, stop_condition):
    """Set simulation parameters."""
    f = h5py.File(filename, 'r+')
    ds_p = f.create_dataset('Parameters', ())
    ds_p.attrs.create('min_cells', data=min_cells)
    ds_p.attrs.create('max_cells', data=max_cells)
    ds_p.attrs.create('max_iters', data=max_iters)
    ds_p.attrs.create('population_size', data=population_size)
    ds_p.attrs.create('spring_force', data=spring_force)
    ds_p.attrs.create('mutation_probability', data=mutation_probability)
    ds_p.attrs.create('stop_condition', data=bytes(stop_condition, TEXT_ENCODING))


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-r', '--residue', 'residues', type=str, multiple=True)
def gene(filename, residues):
    """Define all interesting positions in the gene. Only these will be considered in simulation."""
    f = h5py.File(filename, 'r+')
    infogetter = re.compile(r'([A-Z])([0-9]+)([A-Z]+)')
    positions = []
    options = []
    for r in residues:
        info = infogetter.search(r).groups()
        positions.append(int(info[1]))
        fromopt = info[0]
        toopt = list(info[2])
        toopt.sort()
        options.append([fromopt] + toopt)
    positions, options = zip(*sorted(zip(positions, options)))
    radix = [len(x) for x in options]
    make_rectangular(options, '')
    options = encode_ll(options)
    gp_g = f['Gene']
    gp_g.create_dataset('Radix', data=radix)
    gp_g.create_dataset('Positions', data=positions)
    gp_g.create_dataset('Options', data=options)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('--uni-deathrate', type=float, default = -1.0)
@click.option('--uni-birthrate', type=float, default = -1.0)
@click.option('--uni-minbirthrate', type=float, default = -1.0)
@click.option('-b', '--birthrate', 'birthrates', type=(str, float), multiple=True)
@click.option('-d', '--deathrate', 'deathrates', type=(str, float), multiple=True)
@click.option('-m', '--minbirthrate', 'minbirthrates', type=(str, float), multiple=True)
@click.option('-p', '--wildtype-birthrate', type=float, default = 0.0)
@click.option('-q', '--wildtype-deathrate', type=float, default = 0.0)
@click.option('-w', '--wildtype-minbirthrate', type=float, default = 0.0)
def rate(filename, uni_deathrate, uni_birthrate, uni_minbirthrate, birthrates, deathrates, minbirthrates, wildtype_birthrate, wildtype_deathrate, wildtype_minbirthrate):
    """Define reproduction rates and death rates for cells"""
    f = h5py.File(filename, 'r+')
    gp_rates = f['Rates']
    radix = list(f['Gene/Radix'])
    positions = list(f['Gene/Positions'])
    options = decode_ll(list(f['Gene/Options']))
    # Write unified death/birthrates to file if they are specified.
    # This also server as a marker whether to use them [unified]
    if uni_deathrate >= 0.0:
        gp_rates.attrs.create('unified_death_rate', data=uni_deathrate)
    if uni_birthrate >= 0.0:
        gp_rates.attrs.create('unified_reproduction_rate', data=uni_birthrate)
    if uni_minbirthrate >= 0.0:
        gp_rates.attrs.create('unified_minimum_reproduction_rate', data=uni_minbirthrate)

    genes = []
    birthr = []
    deathr = []
    minbirthr = []

    # FIXME this will not work for variable deathrate
    # TODO implement
    genes.append([0]*len(radix))
    birthr.append(wildtype_birthrate)
    deathr.append(wildtype_deathrate)
    minbirthr.append(wildtype_minbirthrate)

    for r in birthrates:
        g = [0]*len(radix)
        mutant = r[0].split()
        rate = r[1]
        if len(mutant) == 1:
            infogetter = re.compile(r'([A-Z])([0-9]+)([A-Z])')
            info = infogetter.search(mutant[0]).groups()
            position_index = positions.index(int(info[1]))
            option_index = options[position_index].index(info[2])
            g[position_index] = option_index
            # single_mutants[position_index][option_index] = ic50
        else:
            # Set the correct genotype
            # TODO: investigate maybe splitting into a function, as this is reuse:y
            infogetter = re.compile(r'([A-Z])([0-9]+)([A-Z])')
            for m2 in mutant:
                info = infogetter.search(m2).groups()
                position_index = positions.index(int(info[1]))
                option_index = options[position_index].index(info[2])
                g[position_index] = option_index
            # unique_mutants.append((g, ic50))
        genes.append(g)
        birthr.append(rate)

    print(birthrates)
    print(genes)
    print(birthr)

    gp_rates.create_dataset('Genes', data=genes)
    gp_rates.create_dataset('Reproduction', data=birthr)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-n', '--drugname', 'drugname', type=str)
@click.option('-w', '--wildtype', type=float)
@click.option('-m', '--mutant', 'mutants', type=(str, float), multiple=True)
@click.option('-x', '--max_concentration', type=float, default=-1.0)
@click.option('-i', '--min_concentration', type=float, default=0.0)
@click.option('-s', '--ss_concentration', type=float, default=-1.0)
def drug(filename, drugname, wildtype, mutants, max_concentration, min_concentration, ss_concentration):
    """Define a drug. Remember to give an IC50 for at least all SNVs and the wildtype"""
    f = h5py.File(filename, 'r+')
    radix = list(f['Gene/Radix'])
    positions = list(f['Gene/Positions'])
    options = decode_ll(list(f['Gene/Options']))
    unique_mutants = [] # List of anything that is not a single mutant (and has a known ic50)
    # single mutants are listed in a matrix (matches options layout)
    single_mutants = [[0.0]*max(radix) for __ in radix]
    unique_mutants.append(([0]*len(radix), wildtype))
    for m in mutants:
        g = [0]*len(radix)
        mutant = m[0].split()
        ic50 = m[1]
        if len(mutant) == 1:
            infogetter = re.compile(r'([A-Z])([0-9]+)([A-Z])')
            info = infogetter.search(mutant[0]).groups()
            position_index = positions.index(int(info[1]))
            option_index = options[position_index].index(info[2])
            single_mutants[position_index][option_index] = ic50
        else:
            # Set the correct genotype
            # TODO: investigate maybe splitting into a function, as this is reuse:y
            infogetter = re.compile(r'([A-Z])([0-9]+)([A-Z])')
            for m2 in mutant:
                info = infogetter.search(m2).groups()
                position_index = positions.index(int(info[1]))
                option_index = options[position_index].index(info[2])
                g[position_index] = option_index
            unique_mutants.append((g, ic50))
    gp_d = f['Drugs'].create_group(drugname)
    gp_d.create_dataset('Single', data=single_mutants)
    gp_d.create_dataset('Unique_genes', data=[x[0] for x in unique_mutants])
    gp_d.create_dataset('Unique_IC50s', data=[x[1] for x in unique_mutants])
    gp_d.attrs.create('max_concentration', data=max_concentration)
    gp_d.attrs.create('min_concentration', data=min_concentration)
    gp_d.attrs.create('ss_concentration', data=ss_concentration)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-c', '--cell_types', type=(str, int), multiple=True)
def cell(filename, cell_types):
    """Define starting population of system from dna sequences"""
    print(cell_types)
    f = h5py.File(filename, 'r+')
    genomes = [bytes(x[0], TEXT_ENCODING) for x in cell_types]
    counts = [x[1] for x in cell_types]
    gp_c = f['Cells']
    gp_c.create_dataset('Genomes', data=genomes)
    gp_c.create_dataset('Counts', data=counts)


# TODO make sure the rest of the program runs fine (just without fcell functionality) if imports fail
try:
    from Bio.Seq import Seq
    from Bio.Seq import MutableSeq
    from Bio import SeqIO
    pass
except:
    pass
@main.command()
@click.argument('filename', type=click.Path())
@click.argument('fasta', type=click.Path())
@click.option('-c', '--cell_types', type=(str, int), multiple=True)
@click.option('-s', '--selection', type=int, default=-1)
def fcell(filename, fasta, cell_types, selection):
    """Define starting population from mutcodes (like: A123B or for a double mutant: A123B C123D) and a fasta sequence (automatically aligned)"""
    if 'Bio' not in sys.modules:
        print("Failed when importing biopython -> fcell non functional. Use cell instead")
    else:
        f = h5py.File(filename, 'r+')
        positions = f['Gene/Positions'][:]
        residues = [x[0] for x in decode_ll(f['Gene/Options'][:])]
        print(positions)
        print(residues)
        s = SeqIO.read(fasta, 'fasta').seq
        print(len(s), s.translate(), len(s.translate()))
        ref = Seq(''.join(residues))
        print(positions)
        print(ref)
        print('')

        offsets = []

        offset = -len(s)
        while (True):
            cons = Seq('')
            for p in positions:
                cons += s[p*3 + offset : p*3 + offset + 3]
            tcons = cons.translate()
            if (ndiffs(ref, tcons) == 0):
                print(cons)
                print(tcons)
                print('')
                offsets.append((offset, cons))
            offset += 1
            if offset > len(s) - max(positions)*3 - 3:
                break

        print(offsets)
        if len(offsets) > 1:
            if selection < 0:
                print('No unambiguous match, select a match with -s')
            else:
                offsets = offsets[selection]
        dna = offsets[0][1]
        print('\nWildtype DNA')
        print(dna)

        # Now go through the list of mutants, asking the user which one they prefer when multiple changes can create the requested mutation.
        genomes = []
        counts = []
        for c in cell_types:
            mmcm = c[0]
            count = c[1]
            if mmcm == 'WT':
                genomes.append(bytes(str(dna), TEXT_ENCODING))
                counts.append(count)
                continue
            mutcodes = mmcm.split(' ')
            print(mutcodes)
            print(count)
            new_dna = dna.tomutable()
            for mutcode in mutcodes:
                print(mutcode)
                fr, n, m = re.match(r'([A-Z])(\d+)([A-Z])', mutcode).groups()
                n = int(n)
                # print(f, n, m)
                nix = list(positions).index(n)
                assert(residues[nix] == fr)
                # Consider less awful implementation of the possible changes search
                changes = []
                atcg = ['A', 'T', 'C', 'G']
                for i in range(3):
                    for j in atcg:
                        # newcodon = MutableSeq(dna[nix*3 : nix*3 + 3])
                        newcodon = dna[nix*3 : nix*3 + 3].tomutable()
                        newcodon[i] = j
                        # print(newcodon.__repr__())
                        newcodon = newcodon.toseq()
                        if newcodon.translate() == m:
                            changes.append(newcodon[:])
                print(changes)
                change = changes[0]
                if len(changes) > 1:
                    print("More than one change can create this mutant")
                    print('Choose codon manually')
                    for i, ch in enumerate(changes):
                        print(i, ':', ch)
                    choice = click.prompt('Please enter a valid integer', type=int)
                    change = changes[choice]
                elif len(changes) == 0:
                    print("No single change can alter", f, 'to', m)
                    print("Provide codon manually")
                    choice = click.prompt('Please enter a codon', type=str)
                    change = Seq(choice)
                print(change)
                for i in range(3):
                    new_dna[nix*3 + i] = change[i]
            print(mmcm, 'DNA')
            print(new_dna)
            genomes.append(bytes(str(new_dna), TEXT_ENCODING))
            counts.append(count)

        gp_c = f['Cells']
        gp_c.create_dataset('Genomes', data=genomes)
        gp_c.create_dataset('Counts', data=counts)


if __name__ == "__main__":
    main()
