**Instructions**

* The file "script.sh" is to be submitted via slurm for running the job. Please change it as per the requirements of the cluter on which it is being executed.
* To implement a longer runtime case, please modify the parameter "run" in the "input.relaxation" file.
* Finally, open the output file after the run completes to interpret and compare results one of which is "Total wall time"
* Explanation: This is a standard [NVT ensemble](https://www.vasp.at/wiki/index.php/NVT_ensemble) which is a statistical ensemble that is used to study the properties of a material, e.g., velocities and position of particles, potential, interparticle force, acceleration in an iterative manner w.r.t. time, under constant number of particles N, constant Volume V, and constant temperature.

The file `power_setting_benchmarking` runs the LAMMPS script while cycling through an array of power performance options (say `g`) and directs the output to the file `g-output.txt`.
