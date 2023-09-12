%% Homework #2
%Author: Victoria Subritzky Katz
%Date: 09/12/2023



%Params
mu = 10; 
sigma = 2;
amp = 9;
vo = -5; 
N = [5,10,20,40,80,160,1000];           %different number of sampples

plot_flag = 1;                          %plot flag
confidence_int1 = NaN(length(N),2);     %struct for the confidence interval for method #1
confidence_int2 = NaN(length(N),2);     %struct for the confidence interval for method #2
confidence_int3 = NaN(length(N),2);     %struct for the confidence interval for method #3
confidence_int4 = NaN(length(N),2);     %struct for the confidence interval for method #4

for i = 1:length(N)
    % Get samples
    gaus = @(x,mu,sigma,amp,vo)amp*exp(-(((x-mu).^2)/(2*sigma.^2)))+vo;
    x = linspace(-5,25,N(i));
    samples = gaus(x,mu,sigma,amp,vo);

    %samples = normrnd(mu, sigma, N(i), 1);
    sample_mean = sum(samples)/N(i); %calculate sample mean

    %Method #1
    sem = sigma/sqrt(N(i)); 
    confidence_int1(i,:) = [sample_mean-sem*1.96,sample_mean+sem*1.96]; 

    %Method #2
    confidence_int2(i,:) = [sample_mean-sem*2.045,sample_mean+sem*2.045]; 

    %Method #3
    re_samples = normrnd(mu, sigma, N(i), 1);
    re_sample_mean = sum(re_samples)/N(i); %calculate sample mean
    confidence_int3(i,:) = [re_sample_mean-sem*1.96,re_sample_mean+sem*1.96]; 

    %Method 4
    p = 1; %dont know what to set this value as 
    k = 1; %dont know what to set this value as 

    prob_mu_X = NaN(length(samples)); 
    for j = 1:length(samples)
        prob_mu_X(j) = (p/(k*sqrt(2*pi*sigma^2)))*exp((1/2*sigma^2)*(N(i)^2*sigma^2+N(i)*sample_mean^2-2*N(i)*mu*sample_mean^2+N(i)*mu^2));
        
    end

    sem_bayes = sqrt(sum((sample_mean-samples).^2)/N(i)); 
    confidence_int4(i,:) = [sample_mean-sem_bayes*1.96,sample_mean+sem_bayes*1.96]; 


    if plot_flag == 1
        figure
        error_bar = NaN(2,length(x));
        error_bar(1,:) = confidence_int1(i,1);
        error_bar(2,:) = confidence_int1(i,2);

        shadedErrorBar(x,samples,error_bar)

        % labels
        title(sprintf('Gaussian with %d Samples and Method #1 Confidence Intervals', N(i)))
        xlabel('X-axis');
        ylabel('Y-axis');
        %legend('Simulated', 'Theoretical')

        figure
        error_bar = NaN(2,length(x));
        error_bar(1,:) = confidence_int2(i,1);
        error_bar(2,:) = confidence_int2(i,2);

        shadedErrorBar(x,samples,error_bar)

        % labels
        title(sprintf('Gaussian with %d Samples and Method #2 Confidence Intervals', N(i)))
        xlabel('X-axis');
        ylabel('Y-axis');

        figure
        error_bar = NaN(2,length(x));
        error_bar(1,:) = confidence_int3(i,1);
        error_bar(2,:) = confidence_int3(i,2);

        shadedErrorBar(x,samples,error_bar)

        % labels
        title(sprintf('Gaussian with %d Samples and Method #3 Confidence Intervals', N(i)))
        xlabel('X-axis');
        ylabel('Y-axis');

        figure
        error_bar = NaN(2,length(x));
        error_bar(1,:) = confidence_int4(i,1);
        error_bar(2,:) = confidence_int4(i,2);

        shadedErrorBar(x,samples,error_bar)

        % labels
        title(sprintf('Gaussian with %d Samples and Method #4 Confidence Intervals', N(i)))
        xlabel('X-axis');
        ylabel('Y-axis');

    end
end