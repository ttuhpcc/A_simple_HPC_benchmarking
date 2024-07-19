**Nominal Instructions**

* The script to submit to slurm is "script_nocona.sh" which is present within the directory called "RUN" of within the code directory itself. Modify it suitably in case you're running it on another cluster than the "Nocona" partition at TTU HPCC.
* The script "modules_to_loaded.sh" consists of all the modules that need to be loaded. So please ensure you have replacement modules included in the if running on any other system than the "Nocona" partition at TTU HPCC.
* Use "source" command to execute the modules_to_loaded file in case you're running interactively.
* If running using slurm simply use "sbatch" to submit the script "script_nocona.sh"
* Read the "CPU time" from the outout file to estimate performance.
**Important**
  Follow the following steps before running the code:
  * Please download the initial input file to the code from the following [link](https://texastechuniversity-my.sharepoint.com/:u:/g/personal/sagnik_singha_ttu_edu/EYnNOQSCxCJBsvGFH38yOyoBmt9J62aKAU5jM_zpQhvZvA?e=gvXgft).
  * Unzip the downloaded file.
  * Copy or move the extract file to the directory named "RUN" within the code's working directory.
