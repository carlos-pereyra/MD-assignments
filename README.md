# MD-assignments

To get started with these Molecular Dynamics scripts made for UC Davis's MAE 298 special problems, Spring 2019 quarter, Follow these instructions:

### Download

Install these goodies. Anaconda will be used to install a number of python libraries, while ROOT will have to be installed separately for visualization and plotting results.

* [anaconda](https://conda.io/miniconda.html) - Install helper.
* [ROOT](https://github.com/root-project/root) - Plot helper.

#### Conda Installs

Type all these commands in the terminal.

``` 
conda install -c anaconda numpy
```

Optional*
```
conda install -c anaconda numpy
conda install pandas=0.22.0
conda install -c anaconda seaborn
conda install -c conda-forge missingno
conda install -c conda-forge datapackage-py
conda install -c conda-forge matplotlib 

pip install datapackage==1.5

pip install quilt
conda install -c conda-forge pyarrow

conda install -c anaconda python.app 
```
#### Setup ROOT

Make a install directory somewhere on your computer.

``` cd into/install/directory ```

clone ROOT repository

``` git clone https://github.com/root-project/root.git ```

Get inside this new directory which holds all ROOT source code

``` cd /root```

Make a build directory
``` 
mkdir build
cd build
```

Run cmake and make (will work on mac - the two double dots engage cmake in directory one level behind. Meanwhile if I remember correctly 'j8' sets ROOT to work on 8 cores.

``` 
cmake ../root
make -j8
```

### Using the MD scripts in Practice

0. Create conda environment
``` conda create -n py2 python=2.7 ```


1. Activate conda environment
``` source activate py2```

2. Mount ROOT (tell environment where to find ROOT libraries)
``` source /path/to/wherever/this/is/stored/thisroot.sh```

3. run a script

- To write data, verlet (sim1) or langevin algorithms (sim3)
    - ``` python write_*.py ```   
- To plot results
    - ``` python read_*.py ```

4. Deactivate environment
``` source deactivate ```
