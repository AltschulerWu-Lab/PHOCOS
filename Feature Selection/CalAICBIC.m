function [AIC,BIC]=CalAICBIC(X)




dim=size(X);
p=dim(2);
n=dim(1);
% default number of test to get minimun under differnent random centriods
test_num=10;

    % try differnet tests to find minimun disortion under k_temp clusters
    for test_count=2:test_num
        [KL,~,sumd]=kmeans(X,test_count,'replicates',20);
        for i=1:1:max(KL)
        data=X(find(KL==i),:);
        Gamma=cov(data');
        center=repmat(mean(data),size(data,1),1);  % % 
        dis(i)=-1/2*trace((data-center)'*pinv(Gamma+1e-3)*(data-center))-size(data,1)*p/2*log(2*3.14)-size(data,1)/2*log(det(Gamma)+1e-4); %A small number was added to avoid trival solution
        end
        
        LogVal=sum(dis);
        k=test_count*(p);
         AIC(test_count)=-2*LogVal+2*k;
         BIC(test_count)=-2*LogVal+k*log(n);

    end

