function R=CalLayerRank(Obj)
Ap=Obj;
popsize=size(Ap,1); %population size
R=zeros(1,popsize); %initialize layer rank matrix
i=0;

while ~isempty(Ap)
    i=i+1;
    NonD=Nondominated(Ap);
    Paret=[];
    for j=1:length(NonD)
        Paret(j,:)=Ap(NonD(j),:);
    end
    r1=ismember(Obj,Paret,'rows');
    for k=1:length(r1)
        if r1(k)==1
            R(k)=i;
        end
    end
    Objnew=setdiff(Ap,Paret,'rows');
    Ap=Objnew;
end
end
           

 
