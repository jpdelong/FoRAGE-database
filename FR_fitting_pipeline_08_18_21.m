clear; clc;

tic % start matlab benchmark timer

which_type = 2; % specify whether you are fitting type II or type III (now only set up for 2)
num_bootstraps = 200; % how many bootstrapped data sets to use
output_parameters = NaN(num_bootstraps,10); % pre-allocate matrix for parameter output

parfor i = 1:2598 % run through pipeline fitting
    i
    dataset_info = load(['Dataset_',num2str(i),'.mat'],...
            'dim_run','dens2d_units','dens3d_units','eaten_run','std_error_run','sample_size_run','replenished_run','datatype_run','prey_run_2D','prey_run_3D');
    
    % get everything from cm^2 m^2
    if strcmp(dataset_info.dens2d_units(1),'prey per cm2')
        prey_run_2D = dataset_info.prey_run_2D.*1e4;
    elseif strcmp(dataset_info.dens2d_units(1),'prey per m2')
        prey_run_2D = dataset_info.prey_run_2D; % if units are m2 leave it alone
    end
    % get everything from cm^3 m^3
    if strcmp(dataset_info.dens3d_units(1),'prey per cm3')
        prey_run_3D = dataset_info.prey_run_3D.*1e6;
    elseif strcmp(dataset_info.dens3d_units(1),'prey per m3')
        prey_run_3D = dataset_info.prey_run_3D; % if units are m3 leave it alone
    elseif strcmp(dataset_info.dens3d_units(1),'prey per m3') == 0
        prey_run_3D = nan(size(dataset_info.dens3d_units));
    end
            
    if dataset_info.dim_run == 2
        prey_run = prey_run_2D;
    elseif dataset_info.dim_run == 3
        prey_run = prey_run_3D;
    elseif dataset_info.dim_run == 2.5
            if isnan(prey_run_3D(1)) % takes 2D only when there is no 3D
                prey_run = prey_run_2D;
                        prey_run = prey_run.^(2.5/2); % get these to 2.5D
            else
                prey_run = prey_run_3D;
                        prey_run = prey_run.^(2.5/3); % get these to 2.5D
            end 
    end
    
    % check if ANY eaten levels > than prey levels (not allowed)
    prey_scale_test = dataset_info.eaten_run./prey_run; % divide eaten by offered
    if sum(prey_scale_test > 1) > 0
        prey_scale_adjust = 2*max(dataset_info.eaten_run./prey_run); % rescale is twice that maximum factor of eaten / offered
        prey_run = prey_run * prey_scale_adjust; % effectively stretches prey axis so that eaten < offered
    elseif mean(prey_scale_test) < 1e-5 % also, if foraging is too low relative to prey, reduce the prey
        prey_scale_adjust = 1e-2;
        prey_run = prey_run * prey_scale_adjust; % effectively stretches prey axis so that eaten < offered
    else
        prey_scale_adjust = 1;
    end
    
    % pre-allocate vectors
        rsquares = NaN(num_bootstraps,1); % open an empty vector to put in rsquareds
        rss = NaN(num_bootstraps,1); % open an empty vector to put in rss
        AIC = NaN(num_bootstraps,1); % open an empty vector to put in AIC
        parameters = NaN(num_bootstraps,2); % open an empty vector to put in parameters
    
        num_fits = 0;
        q = 1;
        while num_fits < num_bootstraps % keep fitting if you haven't done them all
            %q
            
            % break for data set type - means or raw data
            if strcmp(dataset_info.datatype_run,'Mean')
                % create modeled data set 
                [prey_reshape, eaten_reshape] = model_dataset(dataset_info.std_error_run,dataset_info.sample_size_run,dataset_info.eaten_run,prey_run);
            elseif strcmp(dataset_info.datatype_run,'Raw data')
                % bootstrap data set 
                BS_indices = randsample(length(prey_run),length(prey_run),'true');
                prey_reshape = prey_run(BS_indices);
                eaten_reshape = dataset_info.eaten_run(BS_indices);
            end

            % set up a fit for type I to get starting value for a
            fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,'Upper',Inf,'StartPoint',1);
                ft = fittype('a*x');
                fitobject = fit(prey_reshape,eaten_reshape,ft,fo);
                starting_a = coeffvalues(fitobject);

            % second break in mean data sets - prey replenished or not
            if strcmp(dataset_info.replenished_run,'N')
                if which_type == 2
                    RPE_gof_ft = fittype( 'x-lambertw(a*h*x*exp(-a*(1-h*x)))/(a*h)', 'independent', 'x', 'dependent', 'y' );
                    opts_RPE = fitoptions( RPE_gof_ft );
                        opts_RPE.Lower = [0 0];
                        opts_RPE.StartPoint = [0.1*starting_a 0.01*1/max(eaten_reshape)];
                        opts_RPE.Upper = [Inf Inf];
                        success = 1;
                elseif which_type == 3
                end
                try
                    [fitresult, gof, output] = fit( prey_reshape, eaten_reshape, RPE_gof_ft, opts_RPE);
                catch
                    success = 0;
                end

            elseif strcmp(dataset_info.replenished_run,'Y')
                if which_type == 2
                    HDE_gof_ft = fittype( 'a*x/(1+a*h*x)', 'independent', 'x', 'dependent', 'y' );
                        opts_HDE = fitoptions( HDE_gof_ft );
                        opts_HDE.Lower = [0 0];
                        opts_HDE.StartPoint = [0.01*starting_a 0.01*1/max(eaten_reshape)];
                        opts_HDE.Upper = [Inf Inf];
                        success = 1;
                elseif which_type == 3
                end
                try
                    [fitresult, gof, output] = fit( prey_reshape, eaten_reshape, HDE_gof_ft, opts_HDE);
                catch
                    success = 0;
                end
            end

            if success == 1 
                parameters(q,:) = coeffvalues(fitresult); % output the parameters
                parameters(q,1) = parameters(q,1).*prey_scale_adjust; % grab the a and rescale
                rsquares(q) = gof.rsquare;
                rss(q) = gof.sse;
                AIC(q) = length(prey_reshape)*log(rss(q)/length(prey_reshape))+2*3;
                if rsquares(q) > 0
                    num_fits = num_fits+1;
                    q = q+1;
                end
            end

            % for space clearance rate, take median, adjust units then dimensions
            a_BS = prctile(parameters(:,1),50);
            a_95CI_low = prctile(parameters(:,1),2.5);
            a_95CI_hi = prctile(parameters(:,1),97.5);
            h_BS = prctile(parameters(:,2),50);
            h_95CI_low = prctile(parameters(:,2),2.5);
            h_95CI_hi = prctile(parameters(:,2),97.5);

        end
        
output_parameters(i,:) = [i num_bootstraps mean(rsquares) mean(AIC) a_BS a_95CI_low a_95CI_hi h_BS h_95CI_low h_95CI_hi];

parsave(i,parameters,rsquares,rss,AIC);
    
end

save('output_parameters'); % stores all the output in a mat file
 
et = toc; % grabs the end point for the benchmark counter

save('et');