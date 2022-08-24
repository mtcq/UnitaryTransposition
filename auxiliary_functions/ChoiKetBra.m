%Function transforms a pair of matrices U and V into the ketbra of their Choi vectors, |U>><<V|
%If a single matrix U is used as input, the output is |U>><<U|

%Input: Matrices U and V
%Output: |U>><<V|

%Second input method:
%Input: Matrix U
%Output: |U>><<U|

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

function dKetBraUV = ChoiKetBra(U,V)
if nargin==1
    V=U;
end
dKetU=U(:);
dBraV=V(:)';
dKetBraUV=dKetU*dBraV;
end