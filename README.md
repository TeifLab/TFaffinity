# TFaffinity

This MATLAB code computes protein binding affinities across sequences where the regions are given in a BED file (in this case Example.bed)

## Prerequisites

This code requires bedtools (https://bedtools.readthedocs.io/) and a reference FASTA file. If a different PWM is to be used, then this should be provided (and the filename changed in CalculateProteinAffinity.m)

## Algorithm

The main loop is in TFaffinity_aggregate.m, and performs the following algorithm:

1. Read through each line of the BED file:
    - Extract the sequence from the reference,
    - Compute the binding affinities using the TRAP algorithm (Roider et al., 2007) and the given PWM,
    - Aggregate the affinity and keep track of the number of profiles.
2. Once the file has been cycled through:
    - Return the sum of the profiles,
    - Return the mean of the profile (mean = sum of profiles / number of profiles containing that location).
  
##  Notes

- This does not save individual profiles, but the code is easily tweaked to return single region profiles.
- The code gives affinities at the *left* edge of the motif, so the profiles need shifting by length (motif) /2 to centre the profiles relative to the centre of the regions.
- The code is currently set up to run on a Linux or Unix system due to the use of system calls in MATLAB out to the running shell. These would need to be changed if run on Windows.

### Use on an HPC facility

The code can easily be extended to the embarrassingly parallel computation of large numbers of regions by splitting the regions into more manageable bed files (each of, say, 1000 regions), running the script on each of these bed files as a separate process on a HPC node, and aggregating and averaging the output files generated. 

In this case, to ensure no overwriting of results, the following file names would need to change to incorporate the task identifier given to the process by the HPC facility:
  - Input bed file name,
  - Temporary file name used as workspace,
  - Output file name.

Please consult your local HPC documentation for details of how to do this, as it would be specific to your facility.


## Credits

Written in MATLAB by geejaytee; PWM_affinity.m is a direct translation of an R script in the tRap affinity package developed from Roider et al., 2007.
