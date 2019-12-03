function X=initial(n_pop,n_var)
for i=1:n_var
    X(:,i)=(up(1,i)-lb(1,i))*rand(n_pop,1)+lb(1,i);
end

end