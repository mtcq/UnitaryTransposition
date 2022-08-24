%Function that projects a matrix into the linear space spanned by Choi operators of non-signalling quantum channels

%Input: matrix X, vector DIM with all the local dimensions oredered as [input1 output1 input2 output2 ... ]
%Output: matrix C inside the linear space spanned by Choi operators of NS quantum channels

%Requires: TR.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

function C=ProjNSchannel(X, DIM)
k=max(size(DIM))/2;
C=X;
for i=1:k
    C=C+TR(C,[2*i-1,2*i],DIM)-TR(C,[2*i],DIM);
end

end