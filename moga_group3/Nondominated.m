function NonD=Nondominated(M)

pop=size(M,1);
for i=1:pop
    ni=0;
    for j=1:pop
        if (M(i,1)>M(j,1))&(M(i,2)>M(j,2))
            ni=ni+1;
            L(1,i)=ni;
        else
            L(1,i)=ni;
        end
    end
end
NonD=find(L==0);
end