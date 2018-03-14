#!/usr/bin/python3


import h5py
import click
from matplotlib import pyplot as plt
from cycler import cycler
import re
import statistics
import math
import numpy as np
from itertools import islice
import scipy.stats as stats


TEXT_ENCODING = 'utf-8'


def get_gene(id, radix):
    gene = []
    for x in radix:
        gene.append(id%x)
        id //= x
    return gene

def get_id(gene, radix):
    id = 0
    x = 1
    for i in range(len(gene)):
        id += gene[i] * x
        x *= radix[i]
    return id


def readable_id(id, radix, options, position):
    if id == 0:
        return 'WT'

    gene = get_gene(id, radix)
    if gene == [0]*len(radix):
        return 'WT'
    # print(gene)
    mutcodes = []
    for i, g in enumerate(gene):
        if g != 0:
            mc = options[i][0].decode(TEXT_ENCODING) + str(position[i]) + options[i][g].decode(TEXT_ENCODING)
            mutcodes.append(mc)

    return ' '.join(mutcodes)


# thanks, stackexchange
def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def moving_averages(values, size):
    for selection in window(values, size):
        yield sum(selection) / size


@click.group()
def main():
    pass


@main.command()
@click.argument('filename', type=click.Path())
def penum(filename):
    f = h5py.File(filename, 'r')
    gp_protocol = f['Protocol']
    for i, p in enumerate(gp_protocol):
        print(i, p)


@main.command()
@click.argument('filename', type=click.Path())
@click.argument('id', type=int)
def hid(filename, id):
    """Translate an id to a human readable A123B C456D format."""
    f = h5py.File(filename, 'r')
    gp_gene = f['Gene']
    ds_options = gp_gene['Options']
    ds_position = gp_gene['Positions']
    ds_radix = gp_gene['Radix']

    print(readable_id(id, ds_radix, ds_options, ds_position))


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-i', '--record-id', type=(int, int, int), default = (0, 0, 0))
@click.option('-l', '--type-label-count', type=int, default=5)
@click.option('--liny/--logy', default = True)
@click.option('--nototal/--total', default = False)
def draw(filename, record_id, type_label_count, liny, nototal):
    """Draw a simulation timeline with treatment protocol."""
    f = h5py.File(filename, 'r')
    gp_record = f['Record']
    gp_gene = f['Gene']
    ds_options = gp_gene['Options']
    ds_position = gp_gene['Positions']
    ds_radix = gp_gene['Radix']

    print('Opening:', record_id)
    if (len(gp_record) == 0):
        print('No record in file')
        return
    else:
        gp_this_record = gp_record['/'.join([str(x) for x in record_id])]

        fig, ax = plt.subplots(nrows=2, sharex=True)
        plt.rc('lines', linewidth=2)

        # plot the cell graph
        cindex = gp_this_record['Cells/Index'][:]
        crecord = gp_this_record['Cells/Record'][:]
        xaxis = list(range(len(crecord[0])))

        # we need a list of the cell ids with the highest concentration
        cmaxes = [(id, max(trace)) for id, trace in zip(cindex, crecord)]
        cmaxes = sorted(cmaxes, key=lambda tup: tup[1])
        cmaxes = cmaxes[-type_label_count:]
        cmaxids = [x[0] for x in cmaxes]

        ax[0].set_prop_cycle(cycler('color', ['r', 'g', 'b', 'c', 'm', 'y'])
                             * cycler('linestyle', ['-', '--', ':', '-.']))

        # draw a line representing total population
        if not nototal:
            ax[0].plot(xaxis, [sum(x) for x in zip(*crecord)], label='total', color='lightgrey', linewidth=1)

        for id, trace in zip(cindex, crecord):
            # print labels only for the x ones with highest cell count (at some point)
            if id in cmaxids:
                print(id)
                lbl = readable_id(id, ds_radix, ds_options, ds_position)
                print(lbl)
                ax[0].plot(xaxis, trace, label=lbl)
            else:
                ax[0].plot(xaxis, trace)
        ax[0].set_ylabel('Cell count')
        ax[0].set_xlabel('Timestep')
        ax[0].legend(loc='upper right')
        if not liny:
            ax[0].set_yscale('log')


        # plot the drug dosage graph
        try:
            dindex = gp_this_record['Drugs/Index'][:]
            drecord = gp_this_record['Drugs/Record'][:]

            dcycle = ['r', 'g', 'b', 'c', 'm', 'y']

            for id, trace, col in zip(dindex, drecord, dcycle):
                ax[1].plot(xaxis, trace, label=id.decode(TEXT_ENCODING), color=col)
                ax[1].fill_between(xaxis, 0, trace, alpha=.3, facecolor=col, color = col)
            ax[1].set_ylabel('Drug dose [nM]')
            ax[1].set_xlabel('Timestep')
            ax[1].legend(loc='upper right')
        except:
            pass

        plt.show()


@main.command()
@click.argument('filename', type=click.Path())
def gainfo(filename):
    """Plot genetic algorithm meta-info like fitness~generation."""
    f = h5py.File(filename, 'r')
    gp_record = f['Record']
    n_generations = gp_record.attrs.get('n_generations')[0]
    n_protocols = gp_record.attrs.get('n_protocols')[0]
    n_repeats = gp_record.attrs.get('n_repeats')[0]
    print(n_generations, n_protocols, n_repeats)

    if (len(gp_record) == 0):
        print('No record in file')
        return
    else:
        # Get raw data
        raw_fitness = [[] for __ in range(n_generations)]
        for gn in range(n_generations):
            for pn in range(n_protocols):
                gp_this_record = gp_record['/'.join([str(x) for x in [gn, pn, 0]])]
                fit = gp_this_record.attrs.get('fitness')[0]
                raw_fitness[gn].append(fit)
        # Extract relevant statistics
        # print(raw_fitness)
        fitmean = [statistics.mean(y) for y in raw_fitness]
        fitse = [statistics.stdev(y) / math.sqrt(n_protocols) for y in raw_fitness]
        fitmax = [max(y) for y in raw_fitness]
        fitmin = [min(y) for y in raw_fitness]
        print(fitmean)
        print(fitse)
        print(fitmax)
        print(fitmin)

        x_axis = range(n_generations)
        fig, ax = plt.subplots()
        ax.errorbar(x_axis, fitmean, yerr=fitse, fmt='-o')
        ax.plot(x_axis, fitmax)
        ax.plot(x_axis, fitmin)

        plt.show()


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-g', '--generation', type=int, default=0)
@click.option('-p', '--protocol-count', type=int, default=-1)
@click.option('-a', '--average-mode', type=click.Choice(['mean', 'prob']), default='mean')
def pavg(filename, generation, protocol_count, average_mode):
    """Analyze a set of protocols"""
    f = h5py.File(filename, 'r')
    gp_record = f['Record']
    n_generations = gp_record.attrs.get('n_generations')[0]
    n_protocols = gp_record.attrs.get('n_protocols')[0]
    n_repeats = gp_record.attrs.get('n_repeats')[0]
    if protocol_count < 0:
        protocol_count = n_protocols

    if len(gp_record) == 0:
        print('No record in file')
        return
    else:
        # Get raw data
        raw_fitness = [[] for __ in range(n_generations)]
        for gn in range(n_generations):
            for pn in range(n_protocols):
                gp_this_record = gp_record['/'.join(str(x) for x in [gn, pn, 0])]
                fit = gp_this_record.attrs.get('fitness')[0]
                raw_fitness[gn].append(fit)

        # select protocols to average (add in if elif later)
        protocols = []
        drugs = []
        print('Including protocols:')
        for x in np.array(raw_fitness[generation]).argsort()[-protocol_count:][::-1]:
            print(x, raw_fitness[generation][x])
            gp_this_record = gp_record['/'.join(str(x) for x in [generation, x, 0])]
            dindex = gp_this_record['Drugs/Index'][:]
            drecord = gp_this_record['Drugs/Record'][:]
            dd = {}
            for ix, rc in zip(dindex, drecord):
                dd[ix.decode(TEXT_ENCODING)] = rc[:]
                drugs.append(ix.decode(TEXT_ENCODING))
            protocols.append(dd)
        print(protocols)
        drugs = list(set(drugs))
        print(drugs)

        # average the protocols (add options later)
        avg_protocol = []
        if average_mode == 'mean':
            for d in drugs:
                ap = []
                for i in range(len(protocols[0][d])):
                    ad = []
                    for j in range(len(protocols)):
                        ad.append(protocols[j][d][i])
                    ap.append(statistics.mean(ad))
                avg_protocol.append(ap[:])
        elif average_mode == 'prob':
            for d in drugs:
                ap = []
                for i in range(len(protocols[0][d])):
                    ad = []
                    for j in range(len(protocols)):
                        ad.append(protocols[j][d][i])
                    ad = [1.0 if x > 0.0 else 0.0 for x in ad]
                    ap.append(statistics.mean(ad))
                avg_protocol.append(ap[:])

        print(avg_protocol)

        # draw
        fig, ax = plt.subplots(nrows=2, sharex=True)
        dcycle = ['r', 'g', 'b', 'c', 'm', 'y']

        xaxis = list(range(len(avg_protocol[0])))

        for ix, trace, col in zip(drugs, avg_protocol, dcycle):
            ax[0].plot(xaxis, trace, label=ix, color=col)
            ax[0].fill_between(xaxis, 0, trace, alpha=.3, facecolor=col, color = col)
        ax[0].legend(loc='upper right')

        for ix, trace, col in zip(drugs, avg_protocol, dcycle):
            atr = list(moving_averages(trace, 7))
            atr = [0.0]*3 + atr + [0.0]*3
            print(len(xaxis), len(atr))
            ax[1].plot(xaxis, atr, label=ix, color=col)
        ax[1].legend(loc='upper right')

        plt.show()


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-i', '--protocol-id', type=(int, int), multiple=True)
@click.option('-m', '--mode', type=str)
@click.option('--selected/--all', default = True)
def pcomp(filename, protocol_id, mode, selected):
    """Compare statistic of several protocols"""
    f = h5py.File(filename, 'r')
    gp_record = f['Record']
    # n_generations = gp_record.attrs.get('n_generations')[0]
    # n_protocols = gp_record.attrs.get('n_protocols')[0]
    n_repeats = gp_record.attrs.get('n_repeats')[0]
    max_iterations = f['Parameters'].attrs.get('max_iters')

    if not selected:
        # use all protocols
        print("using all protocols")
        protocol_id = [(0, x) for x in range(len(f['Protocol']))]
    elif len(protocol_id) < 1:
        print("Provide at least one protocol to compare")

    protocol_names = [x for x in f['Protocol']]
    # print(protocol_names)
    # the id's are jumbled, so sort these to match
    protocol_names = [protocol_names[int(x[1])] for x in protocol_id]
    # print(protocol_names)

    if mode == 'wthalf':
        wthalfs = []
        for g, p in protocol_id:
            wtrow = []
            for r in range(n_repeats):
                gp_this_protocol = gp_record['/'.join([str(x) for x in [g, p, r]])]
                wthalf = gp_this_protocol.attrs.get('wt_half')[0]
                wtrow.append(max_iterations if wthalf < 0 else wthalf)
            wthalfs.append(wtrow[:])

        # logwthalfs = [[math.log(x) for x in y if x > 0] for y in wthalfs]
        # squareroot transform is way better. TODO update variable naming
        logwthalfs = [[math.sqrt(x) for x in y if x > 0] for y in wthalfs]
        means = [statistics.mean(x)**2.0 for x in logwthalfs]
        medians = [statistics.median(x) for x in wthalfs]
        print('means\t', means)
        print('medians\t', medians)
        print('medians')
        for a, b in zip(protocol_names, medians):
            print(a, '\t', b)

        print('Repeated t-test comparison of sqrt(wthalf):')
        for a in range(len(protocol_id)):
            for b in range(a, len(protocol_id)):
                if a == b:
                    continue
                try:
                    tt = stats.ttest_ind(logwthalfs[a], logwthalfs[b])
                    print(' ', protocol_names[a], 'to', protocol_names[b], ':', tt, '*' if tt.pvalue < 0.05 else '')
                except:
                    print('  Cannot compare', protocol_names[a], 'to', protocol_names[b])

        fig, ax = plt.subplots(nrows=3, sharex=True)
        # ax[0].boxplot(wthalfs, 0, 'rs', 1)
        # ax[1].boxplot(logwthalfs, 0, 'rs', 1)
        ax[0].violinplot(wthalfs, list(range(len(wthalfs))), points=30, widths=1.0, showmeans=True, showextrema=True, showmedians=True)
        try:
            ax[1].violinplot(logwthalfs, list(range(len(logwthalfs))), points=30, widths=1.0, showmeans=True, showextrema=True, showmedians=True)
        except:
            pass
        ax[2].plot(means)
        ax[2].plot(medians)
        # plt.xticks([y for y in range(len(protocol_id))], protocol_id)
        plt.xticks([y for y in range(len(protocol_id))], protocol_names if protocol_names else protocol_id)
        ax[0].set_xlabel('protocol')
        ax[0].set_ylabel('wthalf')
        ax[1].set_ylabel('sqrt(wthalf)')
        ax[2].set_ylabel('mean(sqrt(wthalf))^2 / median(wthalf)')
        fig.autofmt_xdate()
        plt.show()


if __name__ == "__main__":
    main()
