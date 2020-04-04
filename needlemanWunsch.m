function [score, aln, SC, path, switched] = needlemanWunsch(seq1,seq2, mode, match, mismatch, gap)


len1=length(seq1);
len2=length(seq2);

% switch seqs: image should be rather horizontal than vertical
if len1>len2
    tem = seq1; % store temporarily seq1
    seq1 = seq2;
    seq2 = tem;
    len1=length(seq1);
    len2=length(seq2);
    switched=true;
else
    switched=false;
end


% set default values
if nargin < 3
    mode='similarity';
    match=1;
    mismatch=-1;
    gap=-2;
elseif nargin < 4
    if strcmp(mode,'distance')
        match=-1;
        mismatch=1;
        gap=2;
    else
        match=1;
        mismatch=-1;
        gap=-2;
    end
end

% in case seq or mode would be mixed cases
seq1=upper(seq1);
seq2=upper(seq2);
mode=lower(mode);


% score matrix
SC = zeros(len1+1,len2+1); % preallocate memory


SC(1,2:end) = gap*(1:len2); % populate first row (gap)
SC(2:end,1) = gap*(1:len1); % populate first column (gap)

% traceback matrix
TB = zeros(len1+1,len2+1); % preallocate memory

TB(1,2:end) = 3; % populate first row (gap) up value (index) is 2
TB(2:end,1) = 2; % populate first column (gap) left value (index) is 3


for i=2:len1+1
    for j=2:len2+1
        
        % cases with calculated values
        if seq1(i-1) == seq2(j-1)
            V(1) = SC(i-1,j-1) + match;
        else
            V(1) = SC(i-1,j-1) + mismatch;
        end
        
        % gap cases with calculated values
        V(2) = SC(i-1,j) + gap; % del
        V(3) = SC(i,j-1) + gap; % in
        
        % take max/min and put it in te matrix
        if strcmp(mode,'distance')
            [value,index] = min(V);
        else
            [value,index] = max(V);
        end
        
        SC(i,j) = value;
        TB(i,j) = index;
    end
end

% score is the value in the right bottom corner
score = SC(end,end);

i=size(TB,1); aln1='';
j=size(TB,2); aln2='';

% indexes for displaying the path on graph
path = [1,1]; % starting point

while i>1 || j>1
    path = [[i,j]; path]; %
    switch TB(i,j)
        case 1 % sub/match
            aln1=append(seq1(i-1),aln1);
            aln2=append(seq2(j-1),aln2);
            i = i-1; j = j-1;
        case 2 % indel (up direction)
            aln1=append(seq1(i-1),aln1);
            aln2=append('-',aln2);
            i = i-1;
        case 3 % indel (left direction)
            aln1=append('-',aln1);
            aln2=append(seq2(j-1),aln2);
            j = j-1;
            
    end
 
end

% create char array to visualize alignment
aln(1,:)=aln1;
aln(2,:)=' ';
aln(3,:)=aln2;
aln(2,aln1 == aln2) = '|'; % add the vertical line for each match


end

