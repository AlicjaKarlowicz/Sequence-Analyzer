
%FastaData2 = readFromURL('2H2Z_A','protein');
%FastaData3 = readFromURL('6M03_A','protein');
%[M,switched] = getDotMatrix(FastaData(1).sequence,FastaData(2).sequence,4,1);
%drawDotPlot(M);
%imwrite(1-filM,'dotmat.png');
%[score,aln,sc,path,sw]=needlemanWunsch('A','A');

[score, aln, scMat, path, ind, sw]=smithWaterman('ATAAATCGC','ATACTCGA');
%saveToFastaFormat('test','seq1','seq2',aln,ind);



