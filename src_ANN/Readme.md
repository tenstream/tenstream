### Using Artificial Neural Networks (ANN) to fit LUT
    
Python 2.7x with several modules are needed for the execution.
Required python modules:
* numpy
* netCDF4
* ffnet
* matplotlib (just for plotting routines)

#### Usage:
Just execute Calc_ANN.py with some arguments (use "./Calc_ANN.py --help" for further information):
* --LUT: One or more LUTs to train the ANN
* --ANN_setup: Number of hidden neurons and layers (e.g. "20 20" generates an (number input nodes,20,20,number output nodes)-ANN; number input nodes and number output nodes will be calculated automatically)
* --coeff_type: Use "diffuse", "diff2diff", "dir2diff", or "dir2dir" to calculate the corresponding network. If you use "diffuse" or "diff2diff"/"dir2diff" or "dir2dir" the LUT/s you use for training must be diffuse/direct.
* --test_perc: Percentage of training dataset which is used for testing ( (1.0-test_perc) is used for training). After every training step the dataset is shuffeld, hence the testing-dataset changes after every training.
* --err_inc: If the relative error of the test-dataset after a training step increases more than err_inc compared to the previous training the training will be stopped to prevent overfitting.
* --basename: After every training step the resulting ffnet is saved to "basename_training_step_.net".
* --full (optional): If set a fully connected ANN will be initialized.
* --nproc (optional): Number of processe used for training; Use "ncpu" for maximal possible number; Default is 8.

After every training step the resulting ANN will be saved to "path/to/LUT/". When you start the next time the tenstream-solver set the option "-coeff_mode 1" to use ANNs instead of LUTs.

##### Example:
    ./Calc_ANN.py --LUT path/to/LUT/LUT_8_10.direct.dx67.dy67.pspace.dz20.kabs20.ksca20.g3.phi0.theta0.delta_T_1.000.nc path/to/LUT/LUT_8_10.direct.dx67.dy67.pspace.dz20.kabs20.ksca20.g3.phi0.theta20.delta_T_1.000.nc 
    --ANN_setup 20 20 --coeff_type dir2dir test_perc 0.2 err_inc 0.01 --basename path/to/ANN/ANN_prefix --nproc 4
This line calculates a network:
* using two LUTs for training
* has the shape (6,20,20,64)
* is a dir2dir-ANN
* using 20% of the whole dataset for testing (80% for training)
* stops wenn the error increases more than 1%
* will be saved to "path/to/ANN/ANN_prefix_training_step_.net" after every training step
* using 4 processes for training
        
