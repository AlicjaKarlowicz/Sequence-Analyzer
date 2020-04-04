function [FastaData] = readFromFile(fileName)

fid = fopen(fileName); % get file id
FastaData = parseFasta(fid); % parse fasta fromat

end

