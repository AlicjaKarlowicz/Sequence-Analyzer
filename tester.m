%FastaData = readFromFile('albumin.txt');
FastaData2 = readFromURL('2H2Z_A','protein');
FastaData3 = readFromURL('6M03_A','protein');
%[M,switched] = getDotMatrix(FastaData(1).sequence,FastaData(2).sequence,4,1);
%drawDotPlot(M);
%imwrite(1-filM,'dotmat.png');

[score,aln,sc,path,swi]=needlemanWunsch('A','A');


