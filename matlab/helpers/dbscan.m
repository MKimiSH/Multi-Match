% -------------------------------------------------------------------------
% Function: [class,type]=dbscan(x,k,Eps)
% -------------------------------------------------------------------------
% Input: 
% x - data set (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object 
% Eps - neighborhood radius, if not known avoid this parameter or put []
% -------------------------------------------------------------------------
% Output: 
% class - vector specifying assignment of the i-th object to certain 
% cluster (m,1)
% type - vector specifying type of the i-th object 
% -------------------------------------------------------------------------
% Example of use:
% x=[randn(30,2)*.4;randn(40,2)*.5+ones(40,1)*[4 4]];
% [class,type]=dbscan(x,5,[])
% clusteringfigs('Dbscan',x,[1 2],class,type)
% -------------------------------------------------------------------------

function [class,type]=dbscan(x,k,Eps)

[m,n]=size(x);

if nargin<3 || isempty(Eps)
   [Eps]=epsilon(x,k);   %分析可能的邻域半径Eps，关于函数epsilon，见下面的论述
end

x=[(1:m)' x];  %为每个数据点添加编号，其实是合并两个矩阵
[m,n]=size(x);
type=zeros(1,m);
no=1;
touched=zeros(m,1); %标记某数据点是否已被访问

for i=1:m
    if touched(i)==0;
       ob=x(i,:);
       D=dist(ob(2:n),x(:,2:n));  %使用下面的dist函数生成针对当前点的距离矩阵
       ind=find(D<=Eps);   %见下面的论述

       if length(ind)>1 && length(ind)<k+1
          type(i)=0;
          class(i)=0;
       end
       if length(ind)==1
          type(i)=-1;
          class(i)=-1;  
          touched(i)=1; 
       end

       if length(ind)>=k+1; 
          type(i)=1;
          class(ind)=ones(length(ind),1)*max(no);

          while ~isempty(ind)
                ob=x(ind(1),:);
                touched(ind(1))=1;
                ind(1)=[];
                D=dist(ob(2:n),x(:,2:n));
                i1=find(D<=Eps);

                if length(i1)>1
                   class(i1)=no;
                   if length(i1)>=k+1;
                      type(ob(1))=1;
                   else
                      type(ob(1))=0;
                   end

                   for i=1:length(i1)
                       if touched(i1(i))==0
                          touched(i1(i))=1;
                          ind=[ind i1(i)];   
                          class(i1(i))=no;
                       end                    
                   end
                end
          end
          no=no+1; 
       end
   end
end

i1=find(class==0);
class(i1)=-1;
type(i1)=-1;

%...........................................
function [Eps]=epsilon(x,k)

% Function: [Eps]=epsilon(x,k)
% Input: 
% x - data matrix (m,n); m-objects, n-variables
% k - number of objects in a neighborhood of an object

[m,n]=size(x);

Eps=((prod(max(x)-min(x))*k*gamma(.5*n+1))/(m*sqrt(pi.^n))).^(1/n);

%............................................
function [D]=dist(i,x)

% function: [D]=dist(i,x)                            
% Input: 
% i - an object (1,n)
% x - data matrix (m,n); m-objects, n-variables        
%                                                                 
% Output: 
% D - Euclidean distance (m,1)

[m,n]=size(x);
D=sqrt(sum((((ones(m,1)*i)-x).^2)'));

if n==1
   D=abs((ones(m,1)*i-x))';
end