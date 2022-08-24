%Function that verifies if a symbolic matrix is positive semidefinite using Cholesky decomposition

%Input: a complex symbolic matrix M
%Output: out=1 if PSD, out=0 if not PSD, and a Matrix R such that R'*R==M

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2022

function [out,R] = IsPSDSym(M)
M=sym(M);
if M==M'
    try R=chol(M);
        if R'*R==M
            out=1;
        else
            out=0;
        end
    catch
        out=0;
        R=0;
    end
else
    out=0;
    R=0;
end
end
