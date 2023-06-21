function [t_DSM,Sample_size,score_index,p_all,r_all]=cd_VDI_t_new(brain_id,score_id,cov_id)

Sample_size=size(score_id,1);

    
 k=2;
t_DSM=[];
p_all=[];
r_all=[];

delete(gcp('nocreate'));

for u=1:10000000234
    
try
            parpool(4);
            break;

catch
    u
    disp('cd27: Trying the parpool...');
            continue;
end
end


parfor  i=2:size(brain_id,2)
   
data_lin=[brain_id(:,i),score_id(:,k),cov_id(:,2:end)];
data_lin(any(isnan(data_lin), 2),:) = [];% dele any row NaN

sample_size_lin=size(data_lin,1);
    
    
[r_lin,p_lin]=partialcorr(data_lin(:,1),data_lin(:,2),data_lin(:,3:end));

% r to t

n_lin=sample_size_lin-size(data_lin(:,3:end),2);% sample size-- cov
t_lin=cd_r2t(r_lin,n_lin);


t_DSM=[t_DSM;t_lin];
p_all=[p_all;p_lin];
r_all=[r_all;r_lin];
end
delete(gcp('nocreate'));


score_index=k;

end