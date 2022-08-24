%Function that projects a matrix into the linear space spanned by Choi operators of parallel quantum channels

%Input: matrix X, vector DIM with all the local dimensions oredered as [input1 output1 input2 output2 ... ]
%Output: matrix C inside the linear space spanned by Choi operators of parallel quantum channels

%Requirements: TR.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

function C=ProjParChannel(X, DIM)
k=max(size(DIM))/2;
aux=[];
for i=1:k
    aux=[aux, 2*i];
end
C = X - TR(X,aux,DIM) + trace(X)*eye(prod(DIM))/prod(DIM);
end