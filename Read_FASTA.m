% function for reading in FASTA file
% Length is provided to ensure efficient preallocation of memory
% and should be > length of sequence

function [Seq,seq_length,header] = Read_FASTA(FastaFile,Length)

Seq = zeros(1,Length);
fid4 = fopen(FastaFile,'r');

seq_length = 0;

% If file opens correctly
if fid4>0
    header = fgets(fid4);
    
    % check that it contains a FASTA header
    if header(1)~='>'
        warning(['''' FastaFile ''' is not a FASTA file']);
    else
        header=header(2:end-1);
    
        if uint8(header(end))==13 % Fix for windows-based FASTA files
            % (with CRLF rather than just LF)
            header=header(1:end-1);
        end
    
        disp(['Sequence for ' header ' found']);
        i = 1;
        
        % loop over symbols in sequencce
        while i<=Length
            nextsymbol = fscanf(fid4,'%c',1);
            if ~ischar(nextsymbol) || isempty(nextsymbol)
                warning(['Only ' num2str(i-1) ' symbols found in sequence']);
                break;
            end
            % if FASTA file broken over multiple lines, skip character
            if uint8(nextsymbol)==10 || uint8(nextsymbol)==13  || uint8(nextsymbol)==32
                continue;
            else
                switch lower(nextsymbol)
                    case 'a'
                        Seq(i)=1;
                    case 'c'
                        Seq(i)=2;
                    case 'g'
                        Seq(i)=3;
                    case {'t','u'}
                        Seq(i)=4;
                    case {'x','n'}
                        Seq(i)=5;
                    case 'y'
                        Seq(i)=6;
                    otherwise
                        warning(['Unknown sequence symbol: ' nextsymbol]);
                end
            end
            i=i+1;
        end
 
        seq_length = i-1;
    end
else
    seq_length = 0;
    disp(['''' FastaFile ''' is not readable']);
end

Seq = Seq(1:seq_length);

fclose(fid4);

end
