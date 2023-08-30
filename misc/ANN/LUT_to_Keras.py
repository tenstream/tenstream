#!/usr/bin/env python3
import xarray as xr
import numpy as np
import logging as log

def LUT_to_ANN_input(fname, varname):
    from contextlib import closing
    log.info('Converting LUT input data from {} :: {}'.format(fname, varname))

    with closing(xr.open_dataset(fname)) as D:
        retvars = [ D[d].load().data for d in D.variables if d.endswith(varname) ]
        if len(retvars) != 1:
            raise Exception("Havent found the data I was expecting? :: {}".format(varname))

        dim_names = [k for k in D.keys() if 'values' in k]
        dims = [D[k] for k in dim_names][::-1]

        dim_sizes = [np.size(_) for _ in dims ]

        # Index Space coordinates
        mg = np.mgrid[[slice(_) for _ in dim_sizes]]
        inp_idx = mg.reshape(len(dims),-1).T

        # Physics space coordinates
        grid = np.meshgrid(*dims, indexing='ij')
        inp_phys = np.array(grid).reshape(len(dims),-1).T

        # reverse dims to get it in fortran order, i.e. tau first
        inp_idx = inp_idx[:,::-1]
        inp_phys= inp_phys[:,::-1]

        # convert tau to transmission
        # inp_phys[:,0] = np.exp(-inp_phys[:,0])

        log.info('Converting LUT input data from {} :: {} ... done'.format(fname, varname))
        return inp_idx, inp_phys, retvars[0]


def read_inp_file(inpfile, varname):
    log.info('Reading input data from {} :: {}'.format(inpfile, varname))
    with xr.open_dataset(inpfile) as D:
        data = [ D[_].load().data for _ in ["inp_idx","inp_phys",varname]]
    log.info('Reading input data from {} :: {} ... done'.format(inpfile, varname))
    return data


def average_diffuse_lut(var):
    print("Averaging diffuse coeffs with symmetries")
    v = var.copy().reshape((-1,10,10))
    v[:,0,0] = v[:,1,1] = (v[:,0,0] + v[:,1,1]) / 2
    v[:,0,1] = v[:,1,0] = (v[:,0,1] + v[:,1,0]) / 2
    v[:,0,2] = v[:,0,3] = v[:,0,6] = v[:,0,7] = v[:,1,4] = v[:,1,5] = v[:,1,8] = v[:,1,9] = np.mean(v[:,0,[2,3,6,7]] + v[:,1,[4,5,8,9]], axis=-1)/2
    v[:,1,2] = v[:,1,3] = v[:,1,6] = v[:,1,7] = v[:,0,4] = v[:,0,5] = v[:,0,8] = v[:,0,9] = np.mean(v[:,1,[2,3,6,7]] + v[:,0,[4,5,8,9]], axis=-1)/2
    v[:,2,0] = v[:,3,0] = v[:,4,1] = v[:,5,1] = v[:,6,0] = v[:,7,0] = v[:,8,1] = v[:,9,1] = np.mean(v[:,[2,3,6,7],0] + v[:,[4,5,8,9],1], axis=-1)/2
    v[:,2,1] = v[:,3,1] = v[:,4,0] = v[:,5,0] = v[:,6,1] = v[:,7,1] = v[:,8,0] = v[:,9,0] = np.mean(v[:,[2,3,6,7],1] + v[:,[4,5,8,9],0], axis=-1)/2

    v[:,2,2] = v[:,3,3] = v[:,4,4] = v[:,5,5] = v[:,6,6] = v[:,7,7] = v[:,8,8] = v[:,9,9] = np.mean([ v[:,_,_] for _ in range(2,10)], axis=0)
    v[:,2,3] = v[:,3,2] = v[:,4,5] = v[:,5,4] = v[:,6,7] = v[:,7,6] = v[:,8,9] = v[:,9,8] = np.mean([ v[:,i1,i2] for i1,i2 in zip(range(2,10),(3,2,5,4,7,6,9,8))], axis=0)
    v[:,2,4] = v[:,3,5] = v[:,4,2] = v[:,5,3] = v[:,6,8] = v[:,7,9] = v[:,8,6] = v[:,9,7] = np.mean([ v[:,i1,i2] for i1,i2 in zip(range(2,10),(4,5,2,3,8,9,6,7))], axis=0)

    v[:,2,6] = v[:,2,7] = v[:,3,6] = v[:,3,7] = v[:,4,8] = v[:,4,9] = v[:,5,8] = v[:,5,9] = v[:,6,2] = v[:,6,3] = v[:,7,2] = v[:,7,3] = v[:,8,4] = v[:,8,5] = v[:,9,4] = v[:,9,5] = np.mean([ v[:,i1,i2] for i1,i2 in zip(sorted(list(range(2,10))*2),(6,7,6,7, 8,9,8,9, 2,3,2,3, 4,5,4,5,))], axis=0)
    v[:,2,8] = v[:,2,9] = v[:,3,8] = v[:,3,9] = v[:,4,6] = v[:,4,7] = v[:,5,6] = v[:,5,7] = v[:,6,4] = v[:,6,5] = v[:,7,4] = v[:,7,5] = v[:,8,2] = v[:,8,3] = v[:,9,2] = v[:,9,3] = np.mean([ v[:,i1,i2] for i1,i2 in zip(sorted(list(range(2,10))*2),(8,9,8,9, 6,7,6,7, 4,5,4,5, 2,3,2,3,))], axis=0)
    v[:,2,5] = v[:,3,4] = v[:,4,3] = v[:,5,2] = v[:,6,9] = v[:,7,8] = v[:,8,7] = v[:,9,6] = np.mean([ v[:,i1,i2] for i1,i2 in zip(range(2,10),(5,4,3,2,9,8,7,6))], axis=0)
    return v.reshape((-1,100))


def gen_out_file(inp_idx, inp_phys, varname, var, outfile):
    import xarray as xr
    log.info('Generating output data for {} :: {}'.format(varname, outfile))
    if var.shape[1] == 100:
        var = average_diffuse_lut(var)

    D = xr.Dataset({
        "inp_idx" : xr.DataArray(inp_idx, dims=("Nsample", "Ndim")),
        "inp_phys": xr.DataArray(inp_phys, dims=("Nsample", "Ndim")),
        varname   : xr.DataArray(var, dims=("Nsample", "Ncoeff_{}".format(varname))),
    })
    D.to_netcdf(outfile)

def loss_diff(y_true, y_predict, atol=0, rtol=0):
    import tensorflow as tf
    import tensorflow.math as M
    diff = y_predict - y_true
    lim_atol = M.abs(diff) < atol
    lim_rtol = M.abs(diff/y_true) < rtol
    ret = tf.where(M.logical_and(lim_atol, lim_rtol), 0., diff)
    #if (atol > 0) or (rtol > 0):
    #    tf.print("nonzero diff:",
    #            tf.math.count_nonzero(ret), "/",
    #            tf.size(ret), "(",
    #            (tf.math.count_nonzero(ret)*100)/tf.size(ret, tf.dtypes.int64), "% )")
    return ret

def loss_rmse(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff = loss_diff(y_true, y_predict, **kwargs)
    return M.sqrt( M.reduce_mean( M.square(diff) ) )

def loss_mse(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff = loss_diff(y_true, y_predict, **kwargs)
    return M.reduce_mean( M.square(diff) )

def loss_mae(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff = loss_diff(y_true, y_predict, **kwargs)
    return M.reduce_mean(diff)

def loss_bias(y_true, y_pred):
    import tensorflow as tf
    import tensorflow.math as M

    cc = {9:3, 30:10, 100:10}
    N = cc[y_pred.shape[-1]]
    dst_coeff_y_pred = tf.split(y_pred, axis=1, num_or_size_splits=N)
    dst_coeff_y_true = tf.split(y_true, axis=1, num_or_size_splits=N)

    # then sum over the src dimension
    bias = M.square(M.reduce_sum(dst_coeff_y_pred, axis=0) - M.reduce_sum(dst_coeff_y_true, axis=0))
    mean_bias = M.reduce_mean(bias)/N

    return mean_bias

def custom_loss(y_true, y_pred):
    return loss_mse(y_true, y_pred, atol=1e-0, rtol=1e-2) + loss_bias(y_true, y_pred)

CUSTOM_OBJECTS = {
        "custom_loss" : custom_loss,
        "loss_mae" : loss_mae,
        "loss_mse" : loss_mse,
        "loss_rmse" : loss_rmse,
        "loss_bias" : loss_bias,
        "loss_diff" : loss_diff,
        }

def setup_keras_model(inp, trgt, ident,
        activation='elu',
        activation_output='linear',
        learning_rate=0.01,
        n_hidden=5,
        n_neurons=8,
        optimizer='Adam',
        dropout=None):
    import tensorflow as tf

    name = "{}_f{}_M{}_N{}_drop{}_o{}_opti{}_lr{}".format(
            ident, activation,
            n_hidden, n_neurons,
            dropout, activation_output,
            optimizer, learning_rate)

    log.info('Setup Keras model :: {}'.format(name))

    model = tf.keras.Sequential(name=name)

    act        = getattr(tf.keras.activations, activation)
    output_act = getattr(tf.keras.activations, activation_output)

    def get_N_Neurons(ilay):
        if n_neurons > 0:
            return n_neurons
        elif n_neurons == -1:
            dN_by_dL = float(trgt.shape[1] - inp.shape[1]) / float(n_hidden+1)
            N = inp.shape[1] + np.ceil(dN_by_dL) * (ilay+1)
            print('Linear interpolation of neurons in layer {} ({}) => {}'.format(ilay, dN_by_dL, N))
            return np.int(N)
        elif n_neurons < -1:
            p = np.polynomial.polynomial.Polynomial.fit(
                    [-1,(n_hidden+1)/2,n_hidden],
                    [inp.shape[1], -n_neurons, trgt.shape[1]], 2)
            print('Polynomial interpolation of neurons in layer {} => {}'.format(ilay, p(ilay)))
            return np.int( p(ilay))
        else:
            raise ValueError("Bad n_neurons Values: {}".format(n_neurons))

    kernel_init = "glorot_uniform"

    model.add(tf.keras.layers.Dense(
        get_N_Neurons(0),
        activation=act,
        input_shape=(inp.shape[1],),
        kernel_initializer=kernel_init ))

    if dropout is not None:
        model.add(tf.keras.layers.Dropout(dropout))

    for l in range(1,n_hidden):
        model.add(tf.keras.layers.Dense(
            get_N_Neurons(l),
            activation=act,
            kernel_initializer=kernel_init ))

        if dropout is not None:
            model.add(tf.keras.layers.Dropout(dropout))

    model.add(tf.keras.layers.Dense(
        trgt.shape[1],
        activation=output_act,
        kernel_initializer=kernel_init ))

    opt_type = getattr(tf.keras.optimizers, optimizer)
    opt_args = dict(learning_rate=learning_rate)

    if 'momentum' in inspect.getargspec(opt_type).args:
        opt_args ['momentum'] = .3

    if 'nesterov' in inspect.getargspec(opt_type).args:
        opt_args ['nesterov'] = True

    print (f"Optimizer arguments: {optimizer} {opt_args}")

    opt = opt_type(**opt_args)

    model.compile(loss=custom_loss, optimizer=opt, metrics=[loss_mae, loss_mse, loss_rmse, loss_bias] )

    return model


def train_keras_model(inp, trgt, model, train_split=.1, validation_split=.1, epochs=1000, batch_size=None, outpath=None):
    import numpy as np
    import tensorflow as tf
    import tensorflow_addons as tfa
    import datetime
    import os

    log.info('training Keras model :: {}'.format(model.name))
    callbacks = []

    callbacks.append( tf.keras.callbacks.EarlyStopping(monitor='val_loss', patience=10, restore_best_weights=True) )
    callbacks.append( tfa.callbacks.TimeStopping(seconds=3600*24, verbose=1) )

    print("Check logs with: `tensorboard --logdir logs/ --bind_all`")
    log_dir = os.path.join("logs/", model.name, datetime.datetime.now().strftime("%Y%m%d-%H%M%S"))
    callbacks.append( tf.keras.callbacks.TensorBoard(log_dir=log_dir, histogram_freq=1) )

    callbacks.append( tf.keras.callbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.5, patience=2, min_lr=1e-4) )

    if outpath is not None:
        callbacks.append( tf.keras.callbacks.ModelCheckpoint(filepath=outpath, monitor='val_loss', save_best_only=True) )

    train_mask = np.random.rand(inp.shape[0]) <= train_split
    valid_mask = ~train_mask

    if batch_size is None:
        batch_size = inp.shape[0]

    history = model.fit(
            inp[train_mask],
            trgt[train_mask],
            validation_data=(inp[valid_mask], trgt[valid_mask]),
            epochs=epochs,
            batch_size=batch_size,
            callbacks=callbacks,
            shuffle=True)

    #import timeit
    #run_model = lambda: model.predict(inp)
    #dt = timeit.timeit(run_model, number=1)
    #print("Timing model performance: ", dt)

    return model, history


def tf2fornado(model, phys_input, outpath, accuracy=None, verbose=True):
    import numpy as np
    import xarray as xr

    if model is None:
        D = xr.Dataset()
        D.to_netcdf(outpath)
        return

    layers = [ l for l in model.layers if 'Dropout' not in str(type(l)) ]
    if verbose:
        for i,l in enumerate(layers):
            print("Layer {}: {}, {}".format(i, str(type(l)), l.weights[0].shape))

    D = xr.Dataset()
    D.attrs['Nlayer'] = np.int32(len(layers))
    D.attrs['physical_input'] = np.int32(phys_input)
    D.attrs['keras_name'] = model.name
    D.attrs['accuracy'] = accuracy

    for i, l in enumerate(layers):
        D["w{}".format(i)] = xr.DataArray(l.weights[0].numpy(), dims=("Ninp_{}".format(i), "Nout_{}".format(i)))
        D["b{}".format(i)] = xr.DataArray(l.weights[1].numpy(), dims=("Nout_{}".format(i)))
        D["w{}".format(i)].attrs['activation'] = l.activation.__name__

    D.to_netcdf(outpath)


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
    parser.add_argument('-nt', '--no-train', action="store_true", help='dont run the training cycle')
    parser.add_argument('--train_frac', default=.1, type=float, help='Fraction of LUT which is used to train network')
    parser.add_argument('--learn_rate', default=1e-2, type=float, help='Learning rate for optimizer')
    parser.add_argument('--optimizer', default='Adam', type=str, help='Optimizer Name')
    parser.add_argument('-i', '--load_model', type=str, help='load model from file')
    parser.add_argument('-o', '--save_model', type=str, help='save model to file')
    parser.add_argument('--export', type=str, help='file path to export model for fornado')
    parser.add_argument('-v', '--verbose', action="store_true", help='more verbose output')
    parser.add_argument('--no_gpu', action="store_true", help='dont use the GPU')
    parser.add_argument('-af', '--activation', type=str, help='activation functions in network, e.g. elu, relu, swish')
    parser.add_argument('-ao', '--output_function', type=str, help='activation functions in output layer, e.g. softmax, linear, sigmoid')

    args = parser.parse_args()

    if args.verbose:
        log.basicConfig(level=log.DEBUG)

    if args.no_gpu:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    if args.LUT:
        inp_idx, inp_phys, var = LUT_to_ANN_input(args.LUT, args.varname)
        if not args.dump_input:
            raise ValueError("Have to give dump_input option when reading LUT data")
        gen_out_file(inp_idx, inp_phys, args.varname, var, args.dump_input)
        return

    inp_idx, inp_phys, var = None, None, None
    if args.inp:
        inp_idx, inp_phys, var = read_inp_file(args.inp, args.varname)


    if args.physical_axis:
        log.info("Using physical quantities as input")
        inpvar = inp_phys
    else:
        log.info("Using LUT index sets as input")
        inpvar = inp_idx

    # We write to the file beforehand (bug in xarray if we do it after tf import)
    if args.export:
        tf2fornado(None, None, args.export)

    if args.load_model:
        import tensorflow as tf
        model = tf.keras.models.load_model(args.load_model, custom_objects=CUSTOM_OBJECTS)
        from keras import backend as K
        K.set_value(model.optimizer.lr, args.learn_rate)

    else:
        ident = f"{args.varname}_phys{args.physical_axis}"
        model = setup_keras_model(inpvar, var,
                ident=ident,
                n_hidden=args.Nlayers,
                n_neurons=args.Nneurons,
                optimizer=args.optimizer,
                learning_rate=args.learn_rate,
                activation=args.activation,
                activation_output=args.output_function,
                dropout=args.dropout)


    if not args.no_train:
        train_keras_model(
                inpvar,
                var,
                model,
                train_split=args.train_frac,
                epochs=args.epochs,
                batch_size=args.batch_size,
                outpath=args.save_model)

    # Already done with callback
    #if args.save_model:
    #    model.save(args.save_model, save_format="tf")

    if args.export:
        import tensorflow as tf
        model = tf.keras.models.load_model(args.save_model, custom_objects=CUSTOM_OBJECTS)
        accuracy = model.evaluate(inpvar, var)
        tf2fornado(model, args.physical_axis, args.export, accuracy)


if __name__ == "__main__":
    _main()
