function [FastaData] = readManuallyEntered(seq)
FastaData = struct('id',{},'sequence',{});
FastaData(1).id='MANUALLY_ENTERED';
FastaData(1).sequence=seq;
end

