#!/usr/bin/python3


import h5py
import click
import numpy as np


TEXT_ENCODING = 'utf-8'


@click.group()
def main():
    pass


@main.command()
@click.argument('filename', type=click.Path())
def init(filename):
    """Run first to create a hdf5 file."""
    f = h5py.File(filename, 'w')
    f.create_group('Protocol')
    ds_p = f.create_dataset('Parameters', ())


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-r', '--repeats', type=int)
def setp(filename, repeats):
    """Set general simulation parameters."""
    f = h5py.File(filename, 'r+')
    ds_p = f['Parameters']
    ds_p.attrs.create('repeats', data=repeats)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-p', '--protocols', type=int)
@click.option('-e', '--elites', type=int)
@click.option('-c', '--crossover', type=float)
@click.option('-m', '--mutation', type=float)
@click.option('-r', '--resolution', type=int)
def brep(filename, protocols, elites, crossover, mutation, resolution):
    """Set breeder specific parameters."""
    f = h5py.File(filename, 'r+')
    ds_p = f['Parameters']
    ds_p.attrs.create('protocols', data=protocols)
    ds_p.attrs.create('elites', data=elites)
    ds_p.attrs.create('crossover', data=crossover)
    ds_p.attrs.create('mutation', data=mutation)
    ds_p.attrs.create('breeder_resolution', data=resolution)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-g', '--generations', type=int)
def gamp(filename, generations):
    """Set genetic algorithm specific parameters."""
    f = h5py.File(filename, 'r+')
    ds_p = f['Parameters']
    ds_p.attrs.create('generations', data=generations)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-n', '--name', type=str)
@click.option('-r', '--resolution', type=int, default=1)
@click.option('-l', '--length', type=int)
def newp(filename, name, resolution, length):
    """Init a treatment protocol (for non-managed runs)."""
    f = h5py.File(filename, 'r+')
    gp_protocol = f['Protocol']
    gp_this_protocol = gp_protocol.create_group(name)
    gp_this_protocol.create_dataset('Index', data=[])
    gp_this_protocol.create_dataset('Record', data=[])
    gp_this_protocol.attrs.create('resolution', data=resolution)
    gp_this_protocol.attrs.create('length', data=length)


@main.command()
@click.argument('filename', type=click.Path())
@click.option('-n', '--name', type=str)
@click.option('-d', '--drug', type=str)
@click.option('-q', '--equation', type=str)
@click.argument('params', type=float, nargs=-1)
def edip(filename, name, drug, equation, params):
    """
    options for equation with params:\n
    const: dose\n
    step: 0/1=on/off start, timestep, dose\n
    squarewave: 0/1=on/off start, ontime, offtime, dose\n
    rotate: timesteps
    nullify: start, timesteps
    """
    f = h5py.File(filename, 'r+')
    gp_this_protocol = f['Protocol/' + name]
    ds_index = list(gp_this_protocol['Index'][:])
    ds_record = gp_this_protocol['Record'][:]
    length = gp_this_protocol.attrs.get('length')
    resolution = gp_this_protocol.attrs.get('resolution')

    # If this drug hasn't been included, add it to the index and add a row in the record.
    if bytes(drug, TEXT_ENCODING) not in ds_index:
        ds_index.append(bytes(drug, TEXT_ENCODING))
        del gp_this_protocol['Index']
        gp_this_protocol.create_dataset('Index', data=ds_index)
        if len(ds_record) == 0:
            ds_record = np.array([0.0 for x in range(length)])
            ds_record.shape = (1, length)
        else:
            # ds_record = np.vstack([ds_record, [0.0 for x in range(length)]])
            ds_record = np.vstack([ds_record, np.full((1, length), 0.0, dtype=float)])
    record_row = ds_index.index(bytes(drug, TEXT_ENCODING))

    # this_row = ds_record[record_row]

    # Generate a protocol using the equation
    if equation == 'const':
        constant = params[0]
        ds_record[record_row] += constant
        # ds_record[record_row] = np.array([params[0] + ds_record[record_row][x] for x in range(length)])
    elif equation == 'step':
        constant = params[2]
        steptime = int(params[1])
        if params[0] == 0:
            ds_record[record_row][:steptime] += constant
            # ds_record[record_row] = np.array([params[2] + ds_record[record_row][x] if x < params[1] else ds_record[record_row][x] for x in range(length)])
        elif params[0] == 1:
            ds_record[record_row][steptime:] += constant
            # ds_record[record_row] = np.array([params[2] + ds_record[record_row][x] if x >= params[1] else ds_record[record_row][x] for x in range(length)])
    elif equation == 'squarewave':
        ontime = int(params[1])
        offtime = int(params[2])
        dose = params[3]
        if params[1] + params[2] != 0:
            if params[0] == 0:
                wavearray = np.concatenate([np.full((offtime), 0.0, dtype=float)
                                            , np.full((ontime), dose, dtype=float)])
                # ds_record[record_row] = np.array([params[3] + ds_record[record_row][x] if x % (params[1] + params[2]) < params[1] else ds_record[record_row][x] for x in range(length)])
            elif params[0] == 1:
                wavearray = np.concatenate([np.full((ontime), dose, dtype=float)
                                            , np.full((offtime), 0.0, dtype=float)])
            waveforms = int(length / (ontime + offtime)) + 1
            wavearray = np.tile(wavearray, waveforms)
            wavearray = wavearray[:length]
            ds_record[record_row] += wavearray
                # ds_record[record_row] = np.array([params[3] + ds_record[record_row][x] if x % (params[1] + params[2]) >= params[2] else ds_record[record_row][x] for x in range(length)])
        else:
            print("Warning: squarewave protocol called with 0-duration peaks and valleys -> no effect")
    elif equation == 'rotate':
        rotation = int(params[0])
        np.roll(ds_record[record_row], rotation)
        # ds_record[record_row] = np.array([ds_record[record_row][int((x+params[0]) % length)] for x in range(length)])
    elif equation == 'nullify':
        first = int(params[0])
        last = int(params[0] + params[1])
        ds_record[record_row][first:last] = 0.0
        # ds_record[record_row] = np.array([0.0 if x >= params[0] and x < params[0] + params[1] else ds_record[record_row][x] for x in range(length)])
    else:
        print('Unknown equation:', equation)

    # write modified record to file
    del gp_this_protocol['Record']
    gp_this_protocol.create_dataset('Record', data=ds_record, compression='gzip')


if __name__ == "__main__":
    main()
