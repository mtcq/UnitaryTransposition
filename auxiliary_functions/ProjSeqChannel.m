%Function that projects a matrix into the linear space spanned by Choi operators of sequential quantum channels

%Input: matrix X, vector DIM with all the local dimensions oredered as
%[input1 output1 input2 output2 ... ], vector with order of the parties as
%[1 2 ... k], or [2 1 3], etc. (if no order is given, the order [1 2 ... k] is assumed
%Output: matrix C inside the linear space spanned by Choi operators of sequential quantum channels

%Requires: TR.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 23/08/2022

function C=ProjSeqChannel(X, DIM, order)
nParties=max(size(DIM))/2;
C=X;
if nargin<=2
    aux=[];
    for i=nParties:-1:1
        aux=[aux 2*i];
        C= C - TR(C,aux,DIM) + TR(C,[2*i-1 aux],DIM);
        aux=[aux 2*i-1];
    end
else
    aux=[];
    for i=nParties:-1:1
        j=order(i);
        aux=[aux 2*j];
        C= C - TR(C,aux,DIM) + TR(C,[2*j-1 aux],DIM);
        aux=[aux 2*j-1];
    end
end
end