function population=NewGA(SharedFitness,n_var,n_pop,lb,ub,xp)
SharedFitness=SharedFitness+10;
options=gaoptimset('InitialScores',-SharedFitness,'Generations',2,'Display','off','PopulationSize',n_pop,...
    'CrossoverFraction',0.9,'EliteCount',30,'InitialPopulation',xp);
    [xp,fval,~,output,population]=ga(@fNewGA,n_var,[],[],[],[],lb,ub,[],options);
end

function FNG=fNewGA(x)
FNG=0;
end