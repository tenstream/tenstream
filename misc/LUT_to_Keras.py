#!/usr/bin/env python3
import logging as log

def LUT_to_ANN_input(fname, varname):
    import numpy as np
    import xarray as xr
    from contextlib import closing
    log.info('Converting LUT input data from {} :: {}'.format(fname, varname))

    with closing(xr.open_dataset(fname)) as D:

        dim_names = [k for k in D.keys() if 'values' in k]
        dims = [D[k] for k in dim_names]

        dim_sizes = [np.size(_) for _ in dims[::-1] ]

        # Index Space coordinates
        mg = np.mgrid[[slice(_) for _ in dim_sizes]]
        inp_idx = mg.reshape(len(dims),-1).T

        # Physics space coordinates
        grid = np.meshgrid(*dims[::-1], indexing='ij')
        inp_phys = np.array(grid).reshape(len(dims),-1).T

        retvars = [ D[d].load().data for d in D.variables if d.endswith('S') ]
        if len(retvars) != 1:
            raise Exception("Havent found the data i was expecting?")
        log.info('Converting LUT input data from {} :: {} ... done'.format(fname, varname))
        return inp_idx, inp_phys, retvars[0]


def read_inp_file(inpfile, varname):
    import xarray as xr
    log.info('Reading input data from {} :: {}'.format(inpfile, varname))
    with xr.open_dataset(inpfile) as D:
        return [ D[_].load().data for _ in ["inp_idx","inp_phys",varname]]


def gen_out_file(inp_idx, inp_phys, varname, var, outfile):
    import xarray as xr
    log.info('Generating output data for {} :: {}'.format(varname, outfile))

    D = xr.Dataset({
        "inp_idx" : xr.DataArray(inp_idx, dims=("Nsample", "Ndim")),
        "inp_phys": xr.DataArray(inp_phys, dims=("Nsample", "Ndim")),
        varname   : xr.DataArray(var, dims=("Nsample", "Ncoeff_{}".format(varname))),
    })
    D.to_netcdf(outfile)


def setup_keras_model(inp, trgt, ident,
        activation='elu',
        activation_output='linear',
        learning_rate=0.0001,
        n_hidden=5,
        n_neurons=8,
        optimizer='Adam',
        dropout=None):
    import keras

    name = "{}_f{}_M{}_N{}_drop{}_o{}".format(ident, activation, n_hidden, n_neurons, dropout, activation_output)

    log.info('Setup Keras model :: {}'.format(name))

    model = keras.Sequential(name=name)

    act        = getattr(keras.activations, activation)
    output_act = getattr(keras.activations, activation_output)

    model.add(keras.layers.Dense(n_neurons, activation=act, input_shape=(inp.shape[1],) ))
    if dropout is not None:
        model.add(keras.layers.Dropout(dropout))

    for l in range(1,n_hidden):
        model.add(keras.layers.Dense(n_neurons, activation=act))
        if dropout is not None:
            model.add(keras.layers.Dropout(dropout))

    model.add(keras.layers.Dense( trgt.shape[1], activation=output_act))
    opt = getattr( keras.optimizers, optimizer )(lr=learning_rate)
    model.compile( loss='mse', optimizer=opt, metrics=['mae'] )

    return model


def train_keras_model(inp, trgt, model, train_split=.1, validation_split=.1, epochs=1000, batch_size=None):
    import numpy as np
    import keras
    import datetime
    import os

    log.info('training Keras model :: {}'.format(model.name))

    es_callback = keras.callbacks.EarlyStopping(monitor='val_loss', patience=3)

    print("Check logs with: `tensorboard --logdir logs/ --bind_all`")
    log_dir = os.path.join("logs/", model.name, datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    tensorboard_callback = keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1)

    train_mask = np.random.rand(inp.shape[0]) <= train_split
    valid_mask = ~train_mask

    if batch_size is None:
        batch_size = inp.shape[0]

    history = model.fit(
            inp[train_mask],
            trgt[train_mask],
            validation_split=validation_split,
            validation_data=(inp[valid_mask], trgt[valid_mask]),
            epochs=epochs,
            batch_size=batch_size,
            callbacks=[es_callback, tensorboard_callback])

    return model, history

def _main():
    import argparse
    import os

    def str2bool(v):
        if isinstance(v, bool):
           return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')

    parser = argparse.ArgumentParser(description='Input parameters for Neural Network Generator')
    parser.add_argument('varname', type=str, help='variable which should be trained, e.g. direct.S or direct.T')
    parser.add_argument('-LUT', type=str, help='LUT name which should be converted')
    parser.add_argument('-inp', type=str, help='inp data which has been dumped before')
    parser.add_argument('--physical_axis', const=True, default=False, type=str2bool, nargs='?', help='use physical quantities as input')
    parser.add_argument('-dump_input', type=str, help='output file to dump training and validation data, if set training will not run')
    parser.add_argument('-N', '--Nneurons', default=10, type=int, help='Number of Neurons in each hidden layer')
    parser.add_argument('-M', '--Nlayers', default=2, type=int, help='Number of hidden layers')
    parser.add_argument('-e', '--epochs', default=1000, type=int, help='Number of iterations')
    parser.add_argument('-b', '--batch_size', default=None, type=int, help='Training Batch size, defaults to train dataset size')
    parser.add_argument('-d', '--dropout', default=None, type=float, help='drop out layer strength')
    parser.add_argument('--train_frac', default=.1, type=float, help='Fraction of LUT which is used to train network')
    parser.add_argument('--learn_rate', default=1e-4, type=float, help='Learning rate for optimizer')
    parser.add_argument('-i', '--load_model', type=str, help='load model from file')
    parser.add_argument('-o', '--save_model', type=str, help='save model to file')
    parser.add_argument('-v', '--verbose', action="store_true", help='more verbose output')
    parser.add_argument('--no_gpu', action="store_true", help='dont use the GPU')

    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(level=log.DEBUG)

    if args.no_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    if args.LUT:
        inp_idx, inp_phys, var = LUT_to_ANN_input(args.LUT, args.varname)
    elif args.inp:
        inp_idx, inp_phys, var = read_inp_file(args.inp, args.varname)
    else:
        raise Exception("Have to give either -LUT or -inp option")

    if args.dump_input:
        gen_out_file(inp_idx, inp_phys, args.varname, var, args.dump_input)
        return


    if args.physical_axis:
        log.info("Using physical quantities as input")
        inpvar = inp_phys
    else:
        log.info("Using LUT index sets as input")
        inpvar = inp_idx

    if args.load_model:
        import keras
        model = keras.models.load_model(args.load_model)
    else:
        model = setup_keras_model(inpvar, var, ident=args.varname, n_hidden=args.Nlayers, n_neurons=args.Nneurons, learning_rate=args.learn_rate, dropout=args.dropout)

    train_keras_model(inpvar, var, model, train_split=args.train_frac, epochs=args.epochs, batch_size=args.batch_size)

    if args.save_model:
        model.save(args.save_model)

if __name__ == "__main__":
    _main()
