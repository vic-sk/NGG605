
function laterModelFit_plot(data, fits)

rtsSorted = sort(data);
rts_rec = 1./rtsSorted; 
rts_rec_neg = -rts_rec; 

% Binning for RT and 1/RT plots
rrtBins = 0:0.2:10.0; % cutting off long tail of express saccades

% Compute empirical cumulative RT probabilities
n = length(rtsSorted); 
cumulativeProbabilities = (1:n)./n; 

mu = fits(1);
sigma = fits(2);
%set up the x-axis given the spread of the data
x_axis = linspace(min(rts_rec), max(rts_rec), 100);
%x_axis_neg = linspace(min(rts_rec_neg), max(rts_rec_neg), 100);

%calculations to set up histogram
nbins = 45;
[counts, edges] = histcounts(rts_rec, nbins);
xaxis = edges(1:end-1)+diff(edges);
npdf = counts./trapz(xaxis, counts); %normalize

%plot histogram of 1/RTs
figure
bar(xaxis, npdf);
hold on;
plot(x_axis, normpdf(x_axis, mu, sigma), 'g-', 'LineWidth', 2);
xlabel('1/RTs (sec)');
ylabel('Frequency');



%plot cumulative probabilities
figure
cla reset; hold on;
plot(rts_rec,cumulativeProbabilities,'-ko');
hold on;
plot(fliplr(x_axis),normcdf(x_axis, mu, sigma), 'r-');
xlabel('1/RT (sec)');
ylabel('Cumulative probability');
set ( gca, 'xdir', 'reverse' )


end