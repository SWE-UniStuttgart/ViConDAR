%% Header
%
% fitting power law expoinent to measurements. Given some measurement data
% over heights we can calculate the alpha exponent the fits better the
% values
%
% Vasilis Pettas, Univeersity of Stuttgart, Stuttgart Wind Energy (SWE) 23.05.2019

clc,clear all,close all %#ok<CLALL>

% rng(10); % just fixing random numbers for reproducabilitty

%% Inputs
h_h     = 110;         % hub height meters (reference height)
heights = [70:5:150]; % measurement heights

% since we don't have measurements here examples will be created
speeds = [
    (10.*(heights/h_h).^0.25)+0*randn(1,length(heights));
    (7.*(heights/h_h).^0.14)+0.05*randn(1,length(heights));
    (9.*(heights/h_h).^0.2)+0.4*randn(1,length(heights));
    (15.*(heights/h_h).^0.3)+0.5*randn(1,length(heights));
    (20.*(heights/h_h).^0.3)+0.5*randn(1,length(heights));
    ]; %#ok<*NBRAK>
speeds = repmat(speeds,1,1); % checking scalability

% find uref from heights...
for i = 1:size(speeds,1)      
uref(i,1) = speeds(i,find(heights==h_h)); %#ok<FNDSB>
end

%% just a quick example of how shear profile should look like

% heights_example = [0:5:200];
% h_h_example     = 110;
% alpha_example = 0.2;
% uref_example  = 10;
% example_u = uref_example.*(heights_example/h_h_example).^alpha_example;
% noisyexample_u = example_u+0.7*randn(1,length(example_u));
% awfn_u = awgn(example_u,0.01);
% % awfn_u = example_u+1*randn(1,length(example_u));
% figure, plot (example_u,heights_example)
% hold on
% plot (noisyexample_u,heights_example,'x')
% hold on
% plot (awfn_u,heights_example,'o')
% grid on, xlabel ('Wind speed'), ylabel('Height'), legend ({'analytical solution' 'random noise' 'gausiian noise awgn'},'Location', 'Best')
% xlim ([0 20])
% set(gca,'FontSize',14)
% hold off

%% Fitting procedure with fminsearch 
tic
for i =1:size(speeds,1)
    fcn = @(alphaPL) sum(( uref(i)*(heights/h_h).^(alphaPL) - [speeds(i,:)]).^2); % least square defintion f
    [s,~,~] = fminsearch(fcn, [0.14]);    % Minimise Least-Squares error
    alphaPLfit(i) = s; %#ok<*SAGROW>
end
toc

%% Plotting 
for i =1:size(speeds,1)
    figure
    plot (uref(i).*(heights/h_h).^alphaPLfit(i),heights)
    hold on
    plot (speeds(i,:),heights,'o')  
    grid on, xlabel ('Wind speed'), ylabel('Height'), legend ({'fitted' 'data points'},'Location', 'Best')    
    hold off
    set(gca,'FontSize',14)
end


