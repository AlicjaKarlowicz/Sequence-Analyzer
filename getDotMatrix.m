function [dotMat,switched] = getDotMatrix(seq1,seq2,windowsize,errorlimit)

len1=size(seq1,2);
len2=size(seq2,2);
switched=false;

% switch seqs: image should be rather horizontal than vertical
if len1>len2
    tem = seq1; % store temporarily seq1
    seq1 = seq2;
    seq2 = tem;
    len1=size(seq1,2);
    len2=size(seq2,2);
    switched=true;
end

dotMat = zeros(len1,len2); % pad with zeros

for i = 1:len1 % go through seq1 (rows of dot mat)
    for j =1:len2 % go through seq2 (columns of dot mat)
        if seq1(i) == seq2(j) % compare every nuclobase, if equal
            dotMat(i,j)= 1; % replace zero with 1
        end
    end
end


dotMat = filterDotMat(dotMat,windowsize,errorlimit);

end

