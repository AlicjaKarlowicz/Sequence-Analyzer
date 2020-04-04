function [FastaData] = readFromURL(identifier,seqtype)
% url of ncbi api
url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
data = webread(url,'db',seqtype,'rettype','fasta','id',identifier); % get fasta file by identifier

FastaData = parseFasta(data); % parse fasta format from file from url

end

