      function R=CalLinkPatterns(G,FeatureNum)
cnt=0;
for i=1:1:3
    for j=1:1:3
        num_i=(i-1)*FeatureNum+1;
        num_j=(j-1)*FeatureNum+1;
        mat=[];
        mat=G(num_i:num_i+FeatureNum-1,num_j:num_j+FeatureNum-1)
    R(i,j)=sum(sum(sign(mat)));
    end
end
end