#build crack system
python make_crack_thin_strip.py

#run dynamics
NUM_PROCESSORS=30
mpirun -np $NUM_PROCESSORS python run_dynamics.py

#analyse results
python measure_crack_speed.py