# PX915ReproducibleResult
## Installation instructions (this is *REALLY* important)
The following instructions detail how to install everything you need for this package to work on one of the SCRTP nodes, e.g hetmathsys.

### Loading modules

Once you have logged into an SCRTP node, load up the following modules:

```
module purge; module load GCC/11.2.0 OpenBLAS/0.3.18 CMake/3.22.1 OpenMPI/4.1.1 Python/3.9.6
```

### Making a virtual python environment

The next step is to create a virtual environment, which is essentially makes it much easier to avoid any dependency clashes. Also, it means that once you're done you can just delete the virtual environment and it will get rid of everything.

```
python -m venv myvenv
source myvenv/bin/activate
```

### Installing LAMMPS into venv
1. First, navigate to the directory that this readme is in once you've git cloned this repository, and download LAMMPS locally using git, e.g,
```
git clone -b release https://github.com/lammps/lammps.git mylammps
```
Where `mylammps` is the name of the directory that will be made.

2. There should also be an install_lammps.sh script in this directory. You need to make the script runable with 
```
chmod u+x build_lammps.sh
```

3. Now just run
```
./build_lammps.sh
```
This might take a couple of minutes.

4. Deactivate and re-activate your virtual environment, with e.g,
```
deactivate
source myvenv/bin/activate
```
5. Check the install worked correctly. If you run the following, you should see the output below.
```
python
>>> from lammps import lammps
>>> lmp = lammps()
LAMMPS (15 Sep 2022)
using 1 OpenMP thread(s) per MPI task
```
If you don't see this, then something has gone wrong. Feel free to get in touch with me and I can help you sort it out (installing LAMMPS can be very fiddly). However, if you follow the above exactly step by step then it should work fine. 

### Installing matscipy and it's dependencies
We next need to install matscipy, which contains the bulk of the code that was written over the summer. Navigate to the directory that this readme is in and run the following:
```
git clone https://github.com/libAtoms/matscipy
cd matscipy
pip install .
```
This should install matscipy and it's dependencies.

### Getting the Jupyter Notebook working
You now need to install ipython, and create a kernel for you virtual environment to use the jupyter notebook. Run:
```
pip install ipython=8.5.0
```
Now deactivate and reactivate your virtual environment again, e.g,
```
deactivate
source myvenv/bin/activate
```

Now run
```
ipython kernel install --user --name=myvenv
```
and finally install all the remaining modules
```
pip install nglview h5py 
```

and then (making sure you're in the same directory as this readme), run
```
jupyter-notebook
```
Once the you've opened up the notebook, select 'kernel > change kernel > myvenv'. Now you should be able to run the notebook!