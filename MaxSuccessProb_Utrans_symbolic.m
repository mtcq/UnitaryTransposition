%Script finds the verifies the upper bounds in rigorous way

%Requires: MaxSuccessProb_com_anticom_task.m, MaxSuccessProb_Utrans.m from mtcq

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 19/08/2022

MaxSuccessProb_Utrans;

disp('We now use the computer assisted proof method rigorously certify the following upper bound for causal strategies')

UpperPARsym=sym(879)/sym(1000)
UpperPARfloat=double(UpperPARsym)

UpperSEQsym=sym(936)/sym(1000)
UpperSEQfloat=double(UpperSEQsym)

UpperGENsym=sym(937)/sym(1000)
UpperGENfloat=double(UpperGENsym)

DefineSets_Utrans_symbolic;

CPARsym=sym(CPAR);
CPARsym=ProjParChannel(CPARsym,DIM);
CPARsym=MakePSD(CPARsym);
CPARsym=CPARsym/(trace(CPARsym))*d^2;

if IsPSDSym(UpperPARsym*CPARsym-setPsum/N) && IsPSDSym(UpperPARsym*CPARsym-setMsum/N)
    disp('The upper bound for PAR was certified by a rigorous computer assisted proof method')
else
    disp('The upper bound could not be certified')
end

CSEQsym=sym(CSEQ);
CSEQsym=ProjParChannel(CSEQsym,DIM);
CSEQsym=MakePSD(CSEQsym);
CSEQsym=CSEQsym/(trace(CSEQsym))*d^2;

if IsPSDSym(UpperSEQsym*CSEQsym-setPsum/N) && IsPSDSym(UpperSEQsym*CSEQsym-setMsum/N)
    disp('The upper bound for SEQ was certified by a rigorous computer assisted proof method')
else
    disp('The upper bound could not be certified')
end

CGENsym=sym(CGEN);
CGENsym=ProjParChannel(CGENsym,DIM);
CGENsym=MakePSD(CGENsym);
CGENsym=CGENsym/(trace(CGENsym))*d^2;

if IsPSDSym(UpperGENsym*CGENsym-setPsum/N) && IsPSDSym(UpperGENsym*CGENsym-setMsum/N)
    disp('The upper bound for GEN was certified by a rigorous computer assisted proof method')
else
    disp('The upper bound could not be certified')
end
