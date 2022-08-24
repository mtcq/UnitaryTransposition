%Script defines the sets considered in the discrimination task

%Requires: ChoiKetBra.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

Id=eye(2);
X=[0 1;1 0];
Y=[0 -sqrt(-1);sqrt(-1) 0];
Z=[1 0;0 -1];
XpY=(X+Y)/sqrt(2);
XmY=(X-Y)/sqrt(2);
XpZ=(X+Z)/sqrt(2);
XmZ=(X-Z)/sqrt(2);
YpZ=(Y+Z)/sqrt(2);
YmZ=(Y-Z)/sqrt(2);
ZmY=(Z-Y)/sqrt(2);
ZpY=(Z+Y)/sqrt(2);

%Set channels of setM+ and setM- as Choi opreators
setP(:,:,1)=kron(ChoiKetBra(Id),ChoiKetBra(Id));
setP(:,:,2)=kron(ChoiKetBra(Id),ChoiKetBra(X));
setP(:,:,3)=kron(ChoiKetBra(Id),ChoiKetBra(Z));
setP(:,:,4)=kron(ChoiKetBra(X),ChoiKetBra(Id));
setP(:,:,5)=kron(ChoiKetBra(X),ChoiKetBra(X));
setP(:,:,6)=kron(ChoiKetBra(X),ChoiKetBra(Z));
setP(:,:,7)=kron(ChoiKetBra(XmY),ChoiKetBra(XpY));
setP(:,:,8)=kron(ChoiKetBra(XpY),ChoiKetBra(XmY));
setP(:,:,9)=kron(ChoiKetBra(ZmY),ChoiKetBra(ZpY));
setP(:,:,10)=kron(ChoiKetBra(ZpY),ChoiKetBra(ZmY));

setM(:,:,1)=kron(ChoiKetBra(Y),ChoiKetBra(Id));
setM(:,:,2)=kron(ChoiKetBra(Y),ChoiKetBra(X));
setM(:,:,3)=kron(ChoiKetBra(Y),ChoiKetBra(Z));
setM(:,:,4)=kron(ChoiKetBra(Id),ChoiKetBra(Y));
setM(:,:,5)=kron(ChoiKetBra(X),ChoiKetBra(Y));
setM(:,:,6)=kron(ChoiKetBra(Z),ChoiKetBra(Y));
setM(:,:,7)=kron(ChoiKetBra((Id+sqrt(-1)*Y)/sqrt(2)),ChoiKetBra((Id-sqrt(-1)*Y)/sqrt(2)));
setM(:,:,8)=kron(ChoiKetBra((Id-sqrt(-1)*Y)/sqrt(2)),ChoiKetBra((Id+sqrt(-1)*Y)/sqrt(2)));
