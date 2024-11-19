clear; clc;

to_do = 1:3013;
for i = to_do
    if exist(['DS_',num2str(i),'_fits_type3.mat']) % File exists.
        to_do(to_do==i) = [];
    end
end %1170:1180

%%
tic

num_bootstraps = 200; % how many bootstrapped data sets to use
output_parameters2 = NaN(num_bootstraps,10); % pre-allocate matrix for parameter output
output_parameters3 = NaN(num_bootstraps,13); % pre-allocate matrix for parameter output

parfor i = 1:3013%to_do(8):2237% run through pipeline fitting
%for i = 788%1420:1421%to_do(8):2237%6)%1402:1405%t% run through pipeline fitting
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
    elseif mean(prey_scale_test) < 1e-5
        prey_scale_adjust = 1e-2;
        prey_run = prey_run * prey_scale_adjust; % effectively stretches prey axis so that eaten < offered
    else
        prey_scale_adjust = 1;
    end
    
    % pre-allocate vectors type 2
        rsquares2 = NaN(num_bootstraps,1); % open an empty vector to put in rsquareds
        rss2 = NaN(num_bootstraps,1); % open an empty vector to put in rss
        AIC2 = NaN(num_bootstraps,1); % open an empty vector to put in AIC
        parameters2 = NaN(num_bootstraps,2); % open an empty vector to put in parameters
    
    % pre-allocate vectors type 3
        rsquares3 = NaN(num_bootstraps,1); % open an empty vector to put in rsquareds
        rss3 = NaN(num_bootstraps,1); % open an empty vector to put in rss
        AIC3 = NaN(num_bootstraps,1); % open an empty vector to put in AIC
        parameters3 = NaN(num_bootstraps,3); % open an empty vector to put in parameters
    
    % initiate pipeline starting point
        num_fits = 0;
        q = 1;
        testfitcounter = 1;
        while num_fits < num_bootstraps % keep fitting if you haven't done them all
            q
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
            
            eaten_reshape = eaten_reshape(prey_reshape ~= 0);
            prey_reshape = prey_reshape(prey_reshape ~= 0); 
            eaten_reshape = eaten_reshape(~isnan(eaten_reshape));
            prey_reshape = prey_reshape(~isnan(prey_reshape));          

            % set up a fit for type I to get starting value for a
            fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0,'Upper',Inf,'StartPoint',1);
                ft = fittype('a*x');
                fitobject = fit(prey_reshape,eaten_reshape,ft,fo);
                starting_a = coeffvalues(fitobject);

            % second break - prey replenished or not
            if strcmp(dataset_info.replenished_run,'N')
                % type 2
                    RPE_gof_ft2 = fittype( 'x-lambertw(a*h*x*exp(-a*(1-h*x)))/(a*h)', 'independent', 'x', 'dependent', 'y' );
                        opts_RPE2 = fitoptions( RPE_gof_ft2 );
                        opts_RPE2.Lower = [0 0];
                        opts_RPE2.StartPoint = [starting_a 1/max(eaten_reshape)];
                        opts_RPE2.Upper = [Inf Inf];
                        success = 1;
                try
                    [fitresult2, gof2, output2] = fit( prey_reshape, eaten_reshape, RPE_gof_ft2, opts_RPE2);
                catch
                    success = 0
                end
                
                if success == 1
                    % type 3 ######## a = b*x^z
                        RPE_gof_ft3 = fittype( 'x-lambertw(b*x^z*h*x*exp(-b*x^z*(1-h*x)))/(b*x^z*h)', 'independent', 'x', 'dependent', 'y' );
                            opts_RPE3 = fitoptions( RPE_gof_ft3 );

                            coeffs = coeffvalues(fitresult2); % output the parameters
                            testa = coeffs(1);
                            testh = coeffs(2);

                            opts_RPE3.Lower = [0 0 0];
                            opts_RPE3.StartPoint = [testa testh 0];
                            opts_RPE3.Upper = [Inf Inf 5];
                            success = 1;
                    try
                        if testfitcounter > 1
                            opts_RPE3.Upper = [100*testa 100*testh 6-testfitcounter];
                            if 6-testfitcounter < 1
                                opts_RPE3.Upper = [Inf Inf 0];
                            end
                        end
                        [fitresult3, gof3, output3] = fit( prey_reshape, eaten_reshape, RPE_gof_ft3, opts_RPE3);
                    catch
                        success = 0;
                        testfitcounter = testfitcounter + 1;
                    end
                end

            elseif strcmp(dataset_info.replenished_run,'Y')
                % type 2
                    HDE_gof_ft2 = fittype( 'a*x/(1+a*h*x)', 'independent', 'x', 'dependent', 'y' );
                        opts_HDE2 = fitoptions( HDE_gof_ft2 );
                        opts_HDE2.Lower = [0 0];
                        opts_HDE2.StartPoint = [starting_a 0.1*1/max(eaten_reshape)];
                        opts_HDE2.Upper = [Inf Inf];
                        success = 1;
                try
                    [fitresult2, gof2, output2] = fit( prey_reshape, eaten_reshape, HDE_gof_ft2, opts_HDE2);
                catch
                    success = 0;
                end
                
                if success == 1
                    % type 3
                        HDE_gof_ft3 = fittype( 'b*x^z*x/(1+b*x^z*h*x)', 'independent', 'x', 'dependent', 'y' );
                            opts_HDE3 = fitoptions( HDE_gof_ft3 );

                            coeffs = coeffvalues(fitresult2); % output the parameters
                            testa = coeffs(1);
                            testh = coeffs(2);

                            opts_HDE3.Lower = [0 0 0];
                            opts_HDE3.StartPoint = [testa testh 0];
                            opts_HDE3.Upper = [Inf Inf 3];
                            success = 1;
                    try
                        if testfitcounter > 1
                            opts_RPE3.Upper = [100*testa 100*testh 6-testfitcounter];
                            if 6-testfitcounter < 1
                                opts_RPE3.Upper = [100*testa 100*testh 0];
                            end
                        end
                        [fitresult3, gof3, output3] = fit( prey_reshape, eaten_reshape, HDE_gof_ft3, opts_HDE3);
                    catch
                        success = 0;
                        testfitcounter = testfitcounter + 1;
                    end
                end
            end

            if success == 1
                % type 2
                    parameters2(q,:) = coeffvalues(fitresult2); % output the parameters
                    parameters2(q,1) = parameters2(q,1).*prey_scale_adjust; % grab the a and rescale
                    rsquares2(q) = gof2.rsquare;
                    rss2(q) = gof2.sse;
                    AIC2(q) = length(prey_reshape)*log(rss2(q)/length(prey_reshape))+2*2 + 2*2*(2+1)/(length(prey_reshape)-2-1);
                % type 3
                    parameters3(q,:) = coeffvalues(fitresult3); % output the parameters
                    parameters3(q,1) = parameters3(q,1).*prey_scale_adjust; % grab the a and rescale
                    rsquares3(q) = gof3.rsquare;
                    rss3(q) = gof3.sse;
                    AIC3(q) = length(prey_reshape)*log(rss3(q)/length(prey_reshape))+2*3 + 2*3*(3+1)/(length(prey_reshape)-3-1);
                
                %if rsquares2(q) > 0
                    num_fits = num_fits+1;
                    q = q+1;
                %end
            end

            % for space clearance rate, take median, adjust units then dimensions
            % type 2
                a2_BS = prctile(parameters2(:,1),50);
                a2_95CI_low = prctile(parameters2(:,1),2.5);
                a2_95CI_hi = prctile(parameters2(:,1),97.5);
                h2_BS = prctile(parameters2(:,2),50);
                h2_95CI_low = prctile(parameters2(:,2),2.5);
                h2_95CI_hi = prctile(parameters2(:,2),97.5);
            % type 3
                a3_BS = prctile(parameters3(:,1),50);
                a3_95CI_low = prctile(parameters3(:,1),2.5);
                a3_95CI_hi = prctile(parameters3(:,1),97.5);
                h3_BS = prctile(parameters3(:,2),50);
                h3_95CI_low = prctile(parameters3(:,2),2.5);
                h3_95CI_hi = prctile(parameters3(:,2),97.5);            
                z_BS = prctile(parameters3(:,3),50);
                z_95CI_low = prctile(parameters3(:,3),2.5);
                z_95CI_hi = prctile(parameters3(:,3),97.5);

        end
        
        output_parameters2(i,:) = [i num_bootstraps mean(rsquares2) mean(AIC2) a2_BS a2_95CI_low a2_95CI_hi h2_BS h2_95CI_low h2_95CI_hi];
        output_parameters3(i,:) = [i num_bootstraps mean(rsquares3) mean(AIC3) a3_BS a3_95CI_low a3_95CI_hi h3_BS h3_95CI_low h3_95CI_hi z_BS z_95CI_low z_95CI_hi];

parsave(i,parameters2,rsquares2,rss2,AIC2,2);
parsave(i,parameters3,rsquares3,rss3,AIC3,3);
    
end

%save('output_parameters');  % stores all the output in a mat file
 
et = toc; % grabs the end point for the benchmark counter

save('et');
