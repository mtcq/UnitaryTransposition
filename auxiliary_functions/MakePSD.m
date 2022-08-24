%Function that makes the input matrix M is PSD by adding white noise

%Input: Matrix M
%Output: PSD matrix MPSD, R such that MPSD==R'*R

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2022

function [MPSD,R] = MakePSD(M)
tol=10^(-10);
d=size(M,1);
if isa(M,'sym')
    [PSD, R]=IsPSDSym(M);
    if PSD
        MPSD=M;
    else
        M=(M+M')/2;
        MinEig=min(eig(double(M)));
        MPSD=M+eye(d)*sym(tol+MinEig);
        [PSD, R]=IsPSDSym(MPSD);
        if PSD==0
            error('The output matrix is not PSD!')
        end
    end
else
    M=(M+M')/2;
    try R=chol(M);
        MPSD=M;
    catch
        MinEig=min(eig(M));
        MPSD=M+eye(d)*(tol+MinEig);
        R=chol(MPSD);
    end
end
end
