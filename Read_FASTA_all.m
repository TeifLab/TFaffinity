% function to read in a series of FASTA files given by an input pattern
% if pattern is not provided, then all FASTA files are used
% Note that only the first sequence in each FASTA file is converted
% so files would need to be split into separate FASTA files for all to be read in

function [Seq_DNA] = Read_FASTA_all(varargin)

% if we are given a pattern use that, else read all files in directory
if (nargin==1)
    pattern = varargin{1};
else
    pattern = '*.fa*';
end

FASTA_files = dir(pattern);

if (isempty(FASTA_files))
   warning('No FASTA files found');
   Seq_DNA = [];
   return;
end

FASTA_files = {FASTA_files.name};

% d holds the directories in the path given
d=strsplit(pattern,filesep);

x=length(d);

% d will now hold the containing directory
if ~isempty(d)
    d=strjoin(d(1:(x-1)),filesep);
end

% Seq_DNA is a structure which holds all the relevant data
Seq_DNA = struct('file',FASTA_files,'header',cell(1,length(FASTA_files)),...
    'sequence',cell(1,length(FASTA_files)),...
    'length',cell(1,length(FASTA_files)),...
    'centre_point',cell(1,length(FASTA_files)));

% Cycle through all FASTA files matching pattern
for i=1:length(FASTA_files)
    
    disp(FASTA_files{i});

    % read individual FASTA file
    
    if ~isempty(d)
        [Seq, seq_length, header] = Read_FASTA(strjoin({d,FASTA_files{i}},filesep),1e6);
    else
        [Seq, seq_length, header] = Read_FASTA(FASTA_files{i},1e6);
    end
    
    Seq_DNA(i).header   = header;
    Seq_DNA(i).sequence = Seq;
    Seq_DNA(i).length   = seq_length;
    Seq_DNA(i).centre_point = round(seq_length/2);
   
end

end
