function [] = saveToFile(filename,seq1Name,seq2Name,mode,match,mismatch,gap,score,length,aln)
% save alignment data to proper format

gaps = count(aln(1,:),'-')+count(aln(3,:),'-');
identity = count(aln(2,:),'|');

fid = fopen(append(filename,'.txt'),'wt');
fprintf(fid, '# 1: %s\n# 2: %s\n# Mode: %s\n# Match: %d\n# Mismatch: %d\n# Gap: %d\n# Score: %d\n# Length: %d\n# Identity: %d/%d (%i%%)\n# Gaps: %d/%d (%i%%)\n%s\n%s\n%s',...
    seq1Name,seq2Name,mode, match, mismatch, gap, score, length, identity, length, round(identity*100/length), gaps, length, round(gaps*100/length), aln(1,:), aln(2,:), aln(3,:));
fclose(fid);


end

