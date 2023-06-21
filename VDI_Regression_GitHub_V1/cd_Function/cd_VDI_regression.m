function [h_2_final,intercept,h_2_foundation,intercept_foundation,W]=cd_VDI_regression(VDI,t_square,N,intercept_first)

if nargin == 3
    intercept_first = 1;  
end

% VDI=VDI_MID_SST_Face_T1;
% t_square=t_square;
% N=sample_size;
% intercept_first=1;
%% fix regression 
M=length(VDI);

y=t_square-ones(length(VDI),1)*intercept_first;
X=VDI;



    lin=inv(X'*X);
    
    beta=lin*X'*y;


beta_foundation=beta;
%% h_foundation

h_2_foundation=(beta_foundation*M)/N;
intercept_foundation=intercept_first;
%% second regression source

% W=zeros(M,M);
% for i=1:length(VDI)
%     var_lin=2*(intercept_foundation+beta_foundation*VDI(i))^2;
%     lin=(1/var_lin)*(1/VDI(i));
%     W(i,i)=lin;
% end
% 
% % beta_final
% X=[ones(length(x),1),x];
% lin=inv(X'*W*X);
% beta=lin*X'*W*t_square;
% 
% beta_final=beta(2);
% intercept=beta(1);
%% improve second regression through Q matrix
X=[ones(length(VDI),1),VDI];

W=zeros(M,1);
for i=1:length(VDI)
    var_lin=2*(intercept_foundation+beta_foundation*VDI(i))^2;
    lin=(1/var_lin)*(1/VDI(i));
    W(i)=lin;
end

% Q=X'*W;
Q_lin=zeros(1,M);
for i=1:length(VDI)
   lin=W(i)*VDI(i);
    Q_lin(i)=lin;
end
Q=[W';Q_lin];


lin=inv(Q*X);
beta=lin*Q*t_square;


beta_final=beta(2);
intercept=beta(1);

%% Heritability

h_2_final=(beta_final*M)/N;
end







%%






