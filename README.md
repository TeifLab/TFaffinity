# TFaffinity

This MATLAB code computes protein binding affinities across sequences where the regions are given in a BED file (in this case Example.bed)

## Prerequisites

This code requires bedtools (https://bedtools.readthedocs.io/) and a reference FASTA file.

## Algorithm

1. Read through each line of the BED file
  - Extract the sequence from the reference
  - Compute the binding affinities using TRAP algorithm (Roider et al, 2007)
  - Aggregate the affinity and keep track of the number of profiles
2. Once the file has been cycled through
  - Return the sum of the profiles
  - Return the mean of the profile (mean = sum of profiles / number of profiles containing that location)
  
##  Notes

- This does not save individual profiles, but the code is easily tweaked to return single region profiles.
- The code gives affinities at the *left* edge of the motif, so the profiles need shifting by length (motif) /2 to centre the profiles relative to the centre of the regions.
- If regions are of different lengths, then the centred profiles from the previous note need further centering.
- The code can easily be extended to the embarrassingly parallel (en.wikipedia.org/wiki/Embarrassingly_parallel) computation of large numbers of regions by splitting the regions into more manageable bed files (each of, say, 1000 regions) and aggregating and averaging the output.
