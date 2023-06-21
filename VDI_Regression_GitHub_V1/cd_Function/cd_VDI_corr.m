function [Covariance,intercept,Q_g_foundation,intercept_foundation,W]=cd_VDI_corr(VDI,t_12,N_1,N_2,h_1,h_2,intercept_1,intercept_2,N_s,r)


%% fix regression
M=length(VDI);
intercept_first=(r*N_s)/sqrt(N_1*N_2);

y=t_12-ones(length(VDI),1)*intercept_first;
X=VDI;


    
    lin=inv(X'*X); 
    beta=lin*X'*y;
    
% [beta,bint,r,rint,stats] = regress(y,x);
% beta_foundation=beta;


Q_g_foundation=(beta*M)/sqrt(N_1*N_2);
intercept_foundation=intercept_first;
%% weight for second regression

% coefficient
coefficient_1=(N_1*h_1)/M;
coefficient_2=(N_2*h_2)/M;
Q_g_lin=sum(t_12)/(N_1*N_2);

coefficient_3=(Q_g_lin*sqrt(N_1*N_2))/M;

W=zeros(M,1);
for i=1:length(VDI)
   
var_1=(coefficient_1*VDI(i)+intercept_1)*(coefficient_2*VDI(i)+intercept_2);
var_2=((Q_g_foundation*sqrt(N_1*N_2)*VDI(i))/M)+intercept_foundation;
var_all=var_1+(var_2)^2;

lin=(1/var_all)*(1/VDI(i));    %(1/var)*(1/L_j)

W(i)=lin;
end


%% second regression

x=VDI;
X=[ones(length(x),1),x];


% Q=X'*W;
Q_lin=zeros(1,M);
for i=1:length(VDI)
   lin=W(i)*VDI(i);
    Q_lin(i)=lin;
end
Q=[W';Q_lin];

lin=inv(Q*X);
beta=lin*Q*t_12;


% lin=inv(X'*W*X);
% beta=lin*X'*W*t_12;

beta_final=beta(2);
intercept=beta(1);

Covariance=(beta_final*M)/(sqrt(N_1*N_2));

%% Genetic_Correlation
% r=Covariance/sqrt(h_1*h_2);


end




