function PlotEntropyMap(X,Y)

options = [];

 [COEFF, SCORE] = pca(X);
 X=COEFF(:,1:55);
[a,b]=LDA(Y',options,X);
 
 Projection=X*a(:,1:2);
 color{1}='ws';
 color{2}='wo';
 color{3}='w^';
 color{4}='wd';
 color{5}='wp';
 color{6}='w+';
 color{7}='w+';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The implementation of this function relies on the Statistics Toolbox
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% end

B=mnrfit(Projection,Y);

left_x=min(Projection(:,1));
right_x=max(Projection(:,1));
x_r=right_x-left_x;
up_y=max(Projection(:,2));
down_y=min(Projection(:,2));
y_r=up_y-down_y;


width=200;
length=100;
NewDim1=(Projection(:,1)-left_x)*width/x_r;
NewDim2=(Projection(:,2)-down_y)*length/y_r;

[len,wid]=size(NewDim1);
 Label=zeros(len,max(Y));
 
 for i=1:1:len
 Label(i,Y(i))=1;
 end
 
 
 B=mnrfit([NewDim1,NewDim2],Label);
 
 for j=1:1:width
 for i=1:1:length
 I(i,j)=CalEntropy([j,i],B);
 end
 end
   
 
 imagesc(I)
 
 hold on
 
 
 knum=max(Y);
  for i=1:1:knum
 pos=find(Y==i);
 tmp_x=NewDim1(pos);
 tmp_y=NewDim2(pos);
 
 figure(1)
 plot(tmp_x,tmp_y,color{i},'markersize',25);
 hold on
 end
end


function e=CalEntropy(x,B)

x=[1,x];

val=exp(x*B);
val=[val,1];
p=val./sum(val);

e=sum(-p.*log(p));
end
