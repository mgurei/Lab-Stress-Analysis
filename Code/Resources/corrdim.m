function [D2, Cm, epsilon] = corrdim(x,m,tau,epsilon,sloperange)
% m - embedding dimension
% tau - delay time lag
N=length(x); % length of data vector
k=1:(N-(m-1)*tau); % Number of vector (t)

% Constructed m-dimensional phase space (Embedded dimension) 
xt=x(repmat(k,length(0:tau:(m-1)*tau),1)+repmat((0:tau:(m-1)*tau)',1,length(k)));

% Radii test vector
%epsilon=1:1:500;
% Counter for norm>radii
n=zeros(1,length(epsilon));
for i=1:(N-(m-1)*tau)
    for j=i+1:(N-(m-1)*tau)
       n = n + (epsilon > norm(xt(:,i)-xt(:,j)));
    end
end

Cm=n/((N-(m-1)*tau)*(N-(m-1)*tau-1));
logeps=log2(epsilon);
logCm=log2(Cm);
logCmS=logCm(logCm>sloperange(1) & logCm<sloperange(2));
logepsS=logeps(logCm>sloperange(1) & logCm<sloperange(2));
result=polyfit(logepsS,logCmS,1);
D2=result(1); 
%plot(log2(epsilon),log2(Cm));
end