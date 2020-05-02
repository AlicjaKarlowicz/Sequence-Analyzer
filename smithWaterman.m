function [score, aln, scoreMat, path, indrange, isSwitched] = smithWaterman(seq1,seq2, substituteMat, gap)

len1 = length(seq1);
len2 = length(seq2);


% switch seqs: image should be rather horizontal than vertical
if len1>len2
    tem = seq1; % store temporarily seq1
    seq1 = seq2;
    seq2 = tem;
    len1 = length(seq1);
    len2 = length(seq2);
    isSwitched = true;
else
    isSwitched = false;
end


% set default values
if nargin < 3
    %[AA,AC,AG,AT];[AC,CC,GC,GT];[GA,GC,GG,GT];[TA,TC,TG,TT]
    substituteMat = [[1,-1,-1,-1];[-1,1,-1,-1];[-1,-1,1,-1];[-1,-1,-1,1]];
    gap = -2;
end

% in case seq would be mixed cases
seq1 = upper(seq1);
seq2 = upper(seq2);

% score matrix
scoreMat = zeros(len1+1,len2+1); % preallocate memory

scoreMat(1,2:end) = 0; % populate first row (gap)
scoreMat(2:end,1) = 0; % populate first column (gap)

% traceback matrix
tracebackMat = zeros(len1+1,len2+1); % preallocate memory


tracebackMat(1,2:end) = 4; % populate first row, 4 is the index of 0 case
tracebackMat(2:end,1) = 4; % populate first column


for i = 2:len1+1
    for j = 2:len2+1
        
        % cases with calculated values
        if seq1(i-1) == seq2(j-1)
            
            switch seq1(i-1)
                case 'A'
                    match = substituteMat(1,1);
                case 'C'
                    match = substituteMat(2,2);
                case 'G'
                    match = substituteMat(3,3);
                case 'T'
                    match = substituteMat(4,4);
                otherwise
                    error('Invalid sequence');
            end
            
            v(1) = scoreMat(i-1,j-1) + match;
            
        else % mismatch cases
            
            if seq1(i-1)=='A' || seq2(j-1)=='A'
                
                if seq1(i-1)=='C' || seq2(j-1)=='C'
                    mismatch = substituteMat(1,2);
                    
                elseif seq1(i-1)=='G' || seq2(j-1)=='G'
                    mismatch = substituteMat(1,3);
                    
                elseif seq1(i-1)=='T' || seq2(j-1)=='T'
                    mismatch = substituteMat(1,4);
                    
                end
                
            elseif seq1(i-1)=='C' || seq2(j-1)=='C'
                
                if seq1(i-1)=='G' || seq2(j-1)=='G'
                    mismatch = substituteMat(2,3);
                    
                elseif seq1(i-1)=='T' || seq2(j-1)=='T'
                    mismatch = substituteMat(2,4);
                end
                
            elseif seq1(i-1)=='G' || seq2(j-1)=='G'
                
                if seq1(i-1)=='T' || seq2(j-1)=='T'
                    mismatch = substituteMat(3,4);
                end
            else
                error('Invalid sequence');
            end
            
            v(1) = scoreMat(i-1,j-1) + mismatch;
            
        end
        
        % gap cases with calculated values
        v(2) = scoreMat(i-1,j) + gap; % del
        v(3) = scoreMat(i,j-1) + gap; % in
        
        %zero case
        v(4) = 0;
        
        % take max and put it in te matrix
        [value,index] = max(v);
        
        scoreMat(i,j) = value;
        
        tracebackMat(i,j) = index; % index of V vector determine the direction
        
    end
end

% score is the max value
score = max(scoreMat(:));
[scx,scy] = find(scoreMat == score); % find where are max values, store index

n = length(scx);  % length(scx) = lengths(scy)
aln1 = ''; aln2 = '';


% indexes for displaying the path on graph and for fasta document
path = []; indrange = [];

for n=1:length(scx) % get traceback for every starting point of max value
    
    i = scx(n); j = scy(n);
    indrange = [[i,j]; indrange];
    
    while scoreMat(i,j)~=0 % get path until you meet zero in scoring matrix
        
        path = [[i,j] ;path];
        
        switch tracebackMat(i,j)
            case 1 % sub/match
                aln1 = [seq1(i-1) aln1];
                aln2 = [seq2(j-1) aln2];
                i = i-1; j = j-1;
                
            case 2 % indel (up direction)
                aln1=[seq1(i-1) aln1];
                aln2=['-' aln2];
                i = i-1;
                
            case 3 % indel (left direction)
                aln1 = ['-' aln1];
                aln2 = [seq2(j-1) aln2];
                j = j-1;
                
            case 4 % zero case
                break;
                
        end
    end
    
    indrange = [[i+1,j+1]; [0,0]; indrange];
    aln1 = [' ' aln1];
    aln2 = [' ' aln2];
    
end

% create char array to visualize alignment
aln(1,:) = aln1;
aln(2,:) = ' ';
aln(3,:) = aln2;
aln(2,aln1 == aln2 & aln1 ~= ' ') = '|'; % add the vertical line for each match


end
