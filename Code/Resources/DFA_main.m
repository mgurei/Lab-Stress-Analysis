 
function [Alpha1 Alpha2]=DFA_main(DATA)
% DATA should be a time series of length(DATA) greater than 2000,and of column vector.
%A is the alpha in the paper
%D is the dimension of the time series
%n can be changed to your interest
n=4:1:64;
N1=length(n);
F_n=zeros(N1,1);
 for i=1:N1
     F_n(i)=DFA(DATA,n(i),1);
 end
 plot(log10(n),log10(F_n),'*')
 n=n';
%  plot(log(n),log(F_n));
% xlabel('n')
% ylabel('F(n)')
 %Alpha1
 A=polyfit(log(n(1:12)),log(F_n(1:12)),1);
 Alpha1=A(1);
 %Alpha2
 A=polyfit(log(n(13:end)),log(F_n(13:end)),1);
 Alpha2=A(1);
return