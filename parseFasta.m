function [FastaData] = parseFasta(txt)
% txt is either text from url or file id

% structure with id and seq fields and empty cells
FastaData = struct('id',{},'sequence',{});

ftext = textscan(txt,'%s','delimiter','\n'); % scan text by new line
ftext = ftext{:}; % we get cell array, collect it

% check how many seq are in the file
seqChar = strncmp(ftext,'>',1); % compare first char in ftext with >, returns 1 if true, output is an array
seqNum = sum(seqChar); % number of original sequences in file
seqStarts = [find(seqChar); size(ftext,1)+1]; % get indexes of sequences

% populate FastaData
if seqNum > 0
    
    for n = 1:seqNum
        FastaData(n).id = ftext{seqStarts(n)}(2:end);
        FastaData(n).sequence = strcat(ftext{(seqStarts(n)+1):(seqStarts(n+1)-1)});
    end
    
end

end

