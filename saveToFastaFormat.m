function [] = saveToFastaFormat(filename,seq1Name,seq2Name,aln,ind)
% function writing to file alignment output in given format

fid = fopen(append(filename,'.txt'),'wt');
fprintf(fid, '>%s %s\n%s\n>%s %s\n%s',seq1Name,strrep(strjoin(string(ind(:,1)')),' 0 ','-'),aln(1,:),seq2Name,strrep(strjoin(string(ind(:,2)')),' 0 ','-'),aln(3,:));
fclose(fid);

end

