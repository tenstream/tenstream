from pprint import pprint
import json
import keras as K
#import keras_tuner as Kt
import numpy as np
#np.set_printoptions(threshold=np.inf)
np.set_printoptions(linewidth=np.inf)
import os
import sys
import xarray as xr
import tensorflow as tf
import tensorflow.math as M

def loss_diff(y_true, y_predict):
    import tensorflow as tf
    import tensorflow.math as M
    yp = tf.reshape(y_predict, [-1, Nsrc+1, Ndst])
    #yp = M.maximum(0., M.minimum(1., yp))
    yt = tf.reshape(y_true   , [-1, Nsrc*2+1, Ndst])
    diff = yp[:, :Nsrc, :] - yt[:, :Nsrc, :]
    bias = yp[:, Nsrc, :] - yt[:, Nsrc, :]
    std  = yt[:, Nsrc+1:, :]
    #for i in [0,1,2,3]:
    #    for isrc in range(Nsrc):
    #        tf.print(i, isrc, 'truth:', yt[i,isrc,:], summarize=-1)
    #        tf.print(i, isrc, 'predi:', yp[i,isrc,:], summarize=-1)
    #        tf.print(i, isrc, 'diff :', diff[i,isrc,:], summarize=-1)
    #        tf.print(i, isrc, 'std  :', std[i,isrc,:], summarize=-1)
    #    tf.print(i, 'norm_truth :', yt[i, Nsrc, :], summarize=-1)
    #    tf.print(i, 'norm_predi :', yp[i, Nsrc, :], summarize=-1)
    #    tf.print(i, 'bias :', bias[i,:], summarize=-1)

    return diff, bias, std

def loss_mse(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff, bias, std = loss_diff(y_true, y_predict, **kwargs)
    return M.reduce_mean( M.square(diff) )

def loss_bias(y_true, y_predict, **kwargs):
    diff, bias, std = loss_diff(y_true, y_predict, **kwargs)
    #tf.print('bias variance:', bias_variance, summarize=-1)
    #return M.reduce_mean( M.square(bias) / (2*M.square(tf.reduce_mean(std, axis=-1)) ) )
    return M.reduce_mean( M.abs(bias) )

def loss_rmse(y_true, y_predict, **kwargs):
    import tensorflow as tf
    import tensorflow.math as M
    diff, bias, std = loss_diff(y_true, y_predict, **kwargs)
    #return M.sqrt( M.reduce_mean( M.square(diff) / (2*M.square(std)) ) )
    return M.sqrt( M.reduce_mean( M.square(diff) ))

def lrmse(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff, bias, std = loss_diff(y_true, y_predict, **kwargs)
    return M.reduce_mean( M.square(diff) / (2*M.square(std)), axis=-1)

def lbias(y_true, y_predict, **kwargs):
    import tensorflow.math as M
    diff, bias, std = loss_diff(y_true, y_predict, **kwargs)
    bias_variance = M.reduce_mean( M.square(std), axis=-1 )
    return M.reduce_mean( M.abs(bias) / (2*bias_variance) )

def er(y_true, y_predict, **kwargs):
    return lrmse(y_true, y_predict, **kwargs) / lbias(y_true, y_predict, **kwargs)

def custom_loss(y_true, y_predict, **kwargs):
    import tensorflow as tf
    import tensorflow.math as M
    #tf.print('truth:', y_true.shape, 'predict', y_predict.shape, output_stream=sys.stderr)
    #return loss_rmse(y_true, y_predict) + loss_bias(y_true, y_predict)

    return lrmse(y_true, y_predict, **kwargs) + lbias(y_true, y_predict, **kwargs)

CUSTOM_OBJECTS = {
        "custom_loss" : custom_loss,
        "loss_diff" : loss_diff,
        "loss_mse" : loss_mse,
        "loss_rmse" : loss_rmse,
        "loss_bias" : loss_bias,
        "loss_diff" : loss_diff,
        "lrmse": lrmse,
        "lbias": lbias,
        }

def build_model(inpsize, outsize,
        denselayers=(-64,)*2,
        optimizer="Adam",
        lr=1e-3,
        dropout_inp=0.0,
        dropout_hidden=0.0,
        acti="swish",
        acti_inp="sigmoid",
        acti_out="sigmoid"):
    model = K.models.Sequential()
    model.add(K.Input(shape=(inpsize,)))

    if dropout_inp > 0:
        model.add(K.layers.Dropout(rate=dropout_inp))

    for i, N in enumerate(denselayers):
        Nn = abs(N)
        if N < 0:
            # fit quadratic function to number of inputs, width, outputs
            nlyr = len(denselayers)
            mid = (nlyr-1)//2
            b = np.array([inpsize, Nn, outsize])
            A=np.array([[0,0,1],[mid**2, mid, 1],[nlyr**2, nlyr, 1]])
            qc = np.linalg.solve(A,b)
            f = lambda i, a,b,c: int(a * i**2 + b * i + c)
            Nn = f(i, *qc)

            #N0 = abs(denselayers[0])
            #Nn = int(N0 * (1. - i/nlyr) + outsize)
            print(f"{i=} {N=} {mid=} {nlyr=} -> {Nn=}")

        if i == 0:
            model.add(K.layers.Dense(Nn, activation=acti_inp, use_bias=True))
        else:
            model.add(K.layers.Dense(Nn, activation=acti, use_bias=True))
        if dropout_hidden > 0:
            model.add(K.layers.Dropout(rate=dropout_hidden))
    model.add(K.layers.Dense(outsize, activation=acti_out))
    opt_type = getattr(K.optimizers, optimizer)
    opt = opt_type(learning_rate=lr)
    model.compile(loss=custom_loss, optimizer=opt, metrics=[lrmse, lbias, loss_bias, loss_mse, loss_rmse])
    return model

def tuner_find_best(X, Y, name, verbose=1, force=False):
    project=f"tuner.{name}"

    class MyHyperModel(Kt.HyperModel):
        def build(self, hp):
            model = build_model(
                X.shape[-1],
                Y.shape[-1],
                lr=1e-4, #hp.Float("lr", min_value=1e-4, max_value=1e-2, sampling="log"),
                denselayers=(hp.Choice("lyrwidth", [-6, -20, -30, -50, 6, 20, 32, 128, 512]),) * hp.Choice("nlyr", [1, 2, 3, 4, 6, 8, 12,]),
                optimizer="Adam", #hp.Choice("optimizer", ["Adam", "RMSprop"]),
                dropout_inp=0.0, #hp.Choice("dropout_inp", [0.0, ]),
                dropout_hidden=hp.Choice("dropout_hidden", [0.0, 0.2]),
                acti=hp.Choice("acti", ["relu", "swish"]),
                acti_out=hp.Choice("acti_out", ["linear",]),
            )
            return model

        def fit(self, hp, model, *args, **kwargs):
            return model.fit(
                *args,
                batch_size=hp.Choice("batch_size", [8, 32, 64, 256, 1024, 4096, 32768]),
                **kwargs,
            )

    tunertype = Kt.Hyperband
    tuner = tunertype(
        MyHyperModel(),
        hyperband_iterations=100,
        max_model_size=1e6,
        objective="val_loss",
        max_epochs=100,
        max_consecutive_failed_trials=10,
        directory="./tuner",
        project_name=project,
    )

    my_callbacks = [
        K.callbacks.EarlyStopping(patience=15, verbose=1, restore_best_weights=True),
        # K.callbacks.ModelCheckpoint(filepath=modelname, save_best_only=True, verbose=0),
        # K.callbacks.ReduceLROnPlateau(factor=0.5, patience=10, verbose=1, min_lr=1e-6),
    ]

    modelname = f"best.{project}.h5"

    if not os.path.exists(modelname) or force:
        print(f"Running Tuner search!")
        tuner.search(X, Y, epochs=200, validation_split=0.1, callbacks=my_callbacks)
        best_model = tuner.get_best_models()[0]
        print(f"best model: {best_model}")
        best_model.save(modelname)
    else:
        best_model = K.models.load_model(modelname)

    if verbose:
        best_model.summary()

    def mymodel(X, *args, **kwargs):
        return best_model.predict(X)

    return mymodel

def tf2fornado(model, phys_input, outpath, accuracy=None, verbose=True, inpvars=None):

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
    D.attrs['accuracy'] = json.dumps(accuracy)
    D.attrs['inpvars'] = inpvars

    for i, l in enumerate(layers):
        D["w{}".format(i)] = xr.DataArray(l.weights[0].numpy(), dims=("Ninp_{}".format(i), "Nout_{}".format(i)))
        D["b{}".format(i)] = xr.DataArray(l.weights[1].numpy(), dims=("Nout_{}".format(i)))
        D["w{}".format(i)].attrs['activation'] = l.activation.__name__

    D.to_netcdf(outpath)


var = 'C'

inpdata='../bmcdb/bmc.nc'
inpdata='../bmcdb/bmcdb_ex_ey.annotated.nc'
with xr.open_dataset(inpdata) as D:
    inpvars = ','.join(D['inpvars'].data)
    Nsrc = int(D.attrs[f'Nsrc_{var}'])
    Ndst = int(D.attrs[f'Ndst_{var}'])
    X = D[f'{var}_X'].load().data
    Y = D[f'{var}_Y'].load().data
    X_test = D[f'{var}_X_test'].load().data
    Y_test = D[f'{var}_Y_test'].load().data


class CustomLearningRateScheduler(K.callbacks.Callback):
    def __init__(self, schedule):
        super().__init__()
        self.schedule = schedule

    def on_epoch_begin(self, epoch, logs=None):
        lr = self.model.optimizer.learning_rate
        scheduled_lr = self.schedule(epoch, lr)
        self.model.optimizer.learning_rate = scheduled_lr
        #print(f"\nEpoch {epoch}: Learning rate is {float(np.array(scheduled_lr))}.")



if True:
    import argparse
    parser = argparse.ArgumentParser(description='Input parameters for Neural Network Generator')
    parser.add_argument('-N', '--Nneurons', default=-100, type=int, help='Number of Neurons in each hidden layer')
    parser.add_argument('-M', '--Nlayers', default=4, type=int, help='Number of hidden layers')
    parser.add_argument('-a', '--activation', default='swish', type=str, help='hidden activation functions')
    parser.add_argument('-ai', '--activation_input', default=None, type=str, help='input activation functions')
    parser.add_argument('-ao', '--activation_output', default='sigmoid', type=str, help='output activation functions')
    args = parser.parse_args()

    width, depth, acti, acti_inp, acti_out = args.Nlayers, args.Nneurons, args.activation, args.activation_input, args.activation_output

    if acti_inp is None:
        acti_inp = acti

    outname = f'./S.{acti}.{acti_inp}.{acti_out}.{depth}.{width}.nc'

    print(f"Model path {outname}")
    t = int(X.shape[0]*.75)
    X_train, Y_train, X_val, Y_val = X[:t], Y[:t], X[t:], Y[t:]

    print(f"Shape X_train {X_train.shape}")
    print(f"Shape Y_train {Y_train.shape}")
    print(f"Shape X_val   {X_val.shape}")
    print(f"Shape Y_val   {Y_val.shape}")

    print(f"Input Data size: {X.shape} -> {Y.shape}, Y_predict: {(Nsrc+1)*Ndst}")

    if os.path.exists(f'{outname}.keras'):
        print(f"Loading existing model from {outname}.keras")
        model = K.models.load_model(f'{outname}.keras', custom_objects=CUSTOM_OBJECTS)
        initial_value_threshold = M.reduce_mean(custom_loss(Y_val, model.predict(X_val, batch_size=int(1e5)))).numpy()
        print(f"{initial_value_threshold=}")
    else:
        model = build_model(X_train.shape[-1], (Nsrc+1)*Ndst, denselayers=(depth,)*width, lr=1., acti=acti, acti_inp=acti_inp, acti_out=acti_out)
        initial_value_threshold = None

    Nparam = model.count_params()
    initial_lr = 1e-4 #Nparam * 1e-8 / (Nsrc*Ndst)
    print(f"Have {Nparam} parameters -> lr {initial_lr}")

    def lr_schedule(epoch, lr):
        if epoch < 1:
            return initial_lr
        if (epoch % 50) == 0:
            new_lr = 10 * lr
            print(f"increasing learning rate from {lr} to {new_lr}")
            return new_lr
        return lr

    my_callbacks = [
        K.callbacks.EarlyStopping(patience=20, verbose=1, restore_best_weights=True),
        K.callbacks.ReduceLROnPlateau(factor=0.5, patience=3, verbose=1, min_lr=1e-7),
        CustomLearningRateScheduler(lr_schedule),
        K.callbacks.ModelCheckpoint(f'{outname}.keras', monitor='val_loss', save_best_only=True, verbose=True, initial_value_threshold=initial_value_threshold)
    ]
    fit = dict(
            callbacks=my_callbacks,
            validation_data=(X_val, Y_val),
            epochs=10000,
            shuffle=True,
            )


    model.summary()

    model.fit(X_train, Y_train, batch_size=int(128), **fit)

    for jl in range(3):
        for bs in (128,):
            print(f"Iteration {jl} :: batchsize {bs}")
            #K.backend.set_value(model.optimizer.learning_rate, lr)
            #model.fit(X_train, Y_train, batch_size=int(bs), **fit)

    k = int(X.shape[0]) -5
    print("X         ",  X[k:k+1])
    print("Y (target)\n",  Y[k:k+1].reshape(-1, 21,10))

    Yp = model.predict(X[k:k+1]).reshape(-1, 11,10)
    print("Yp        \n", Yp)

    k = int(X_test.shape[0]) -5
    print("X_test         ", X_test   [k:k+1])
    print("Y_test (target)\n", Y_test   [k:k+1].reshape(-1, 21,10))
    Y_predict = model.predict(X_test[k:k+1]).reshape(-1, 11,10)
    print("Y_predict      \n", Y_predict)#[k])

    accuracy = list(zip(model.metrics_names, model.evaluate(X_test, Y_test)))
    tf2fornado(model, True, outname, accuracy, inpvars=inpvars)

else:
    tuner_find_best(X, Y, f'C')
