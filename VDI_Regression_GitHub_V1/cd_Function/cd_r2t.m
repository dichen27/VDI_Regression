function  t = cd_r2t(r,n)  
% source
% https://www.polyu.edu.hk/mm/effectsizefaqs/effect_size_equations2.html
df=n-2;
t=zeros(length(r),1);
for i=1:length(r)
    
t(i)=sqrt(df)*(r(i)/sqrt(1-r(i)^2));

end

end

