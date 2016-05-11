function Pos=findNearResposeCurve(Latent,DataMat,num) 

for i=1:1:size(Latent,1)
    curvec=Latent(i,:);
    Mat=repmat(curvec,size(DataMat,1),1);
    Dis=sum((Mat-DataMat).*(Mat-DataMat),2);
    
    [~,order]=sort(Dis,'ascend');
    
    Pos(i,1:num)=order(1:num);
    
end