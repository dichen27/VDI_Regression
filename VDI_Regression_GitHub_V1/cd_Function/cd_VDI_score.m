function VDI=cd_VDI_score(data)


data_R_matrix=corr(data,data,'rows','complete');% correlation without NaN
M=size(data_R_matrix,1);
base_0=zeros(1,M);
VDI=zeros(M,1);
  for i=1:M   
%     i
   VDI_lin=pdist2(data_R_matrix(i,:),base_0)^2-1;
   VDI(i)=VDI_lin;
  end
end