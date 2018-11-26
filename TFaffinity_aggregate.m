% Script for parsing a single bed file and generating aggregated profiles
% Note that this does not save individual profiles as it only outputs the 
% sum of the profiles and the mean profile 

% For Windows systems, the lines using "system" command would need to change

% Open bed file - change the file name for each file that you wish to parse
fid8 = fopen(['Example.bed']);

% Set up temporary file for processing later
temp_file = ['temp-file'];

% Change this for largest possible sequence in the BED file
max_length = 4000;

i = zeros(1,max_length);

running_total = zeros(1,max_length);


% Set multiplier for TRAP algorithm
global K0
K0 = 1e9;

% Read through the file if it exists
if (fid8>-1)
   
   % While still lines in BED file
   while ~feof(fid8)
      
      % read in line from BED file
      line = fgets(fid8);

      % extract single line
      system (['echo "' line '" > ' temp_file '.bed']);

      % Use BedTools to obtain sequence for single line 
      % Change the .fa file for the reference FASTA file
      system (['bedtools getfasta -fi mm9.fa -bed ' temp_file '.bed -fo ' temp_file '.fa']);
  
      % Read in fasta file to sequence format
      Seq_DNA=Read_FASTA_all([temp_file '.fa']);

      disp(Seq_DNA.header);

      profiles = zeros(1,length(Seq_DNA.sequence));

      % Compute protein affinity profile
      Protein = ComputeProteinAffinity(Seq_DNA.sequence);

      profiles(1:length(Protein)) = Protein;

      % Compute running total and running number of profiles per location
      running_total(1:length(Protein)) = running_total(1:length(Protein)) + profiles(1:length(Protein));

      i(1:length(Protein)) = i(1:length(Protein)+1;
   end

   % Compute averages - note that some regions may be shorter, in which case
   % i(j) returns the number of profiles starting at that position j

   average = running_total ./ i;

end

% Output averages:
% note that the output file consists of affinity at the *left* edge of the
% binding site, so need shifting by half the length
% to centre the profile on the centre of the binding site

save(['average.mat'],'i','average','running_total');

% Clean up temporary files
system(['rm ' temp_file '.bed ' temp_file '.fa' ]);




