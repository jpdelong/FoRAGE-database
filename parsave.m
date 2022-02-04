function parsave(i,parameters,rsquares,rss,AIC)

save(['DS_',num2str(i),'_type2.mat'],...
        'parameters','rsquares','rss','AIC');
    
end