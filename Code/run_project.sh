#/bin/bash
#Author: Ruth Keane (ruth.keane19@imperial.ac.uk)
#Script: run_Project.sh
#Desc: runs project files to simulate and plot
#Date: August 2020
#from Code directory
mkdir ../Results
mkdir ../Results/fromfunction
mkdir ../Results/fromfunction/cors
mkdir ../Writeup 

echo "Running Analytical Sensitivity Analysis "
python3 datamanipulation/pythonsensitivity.py
echo "Run main integration -1 and 10 simulations"
cd Model
Rscript integration_Vectorisation_cov_var.R
echo "Run Plotting and Analysis"
cd ../writeup
Rscript PlotsFinal.R
echo "Finished!"
### note , not tested yet

# to do this - copy across most recent version of all files. # link to repo and use gitignore to lose data_mod, results and writeup. #check runs
#git ignore- nothing from results 

#also want to move all data modifications to a folder called data_mod (change in files that use and put that in gitignore)