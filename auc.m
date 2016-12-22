function [ac,R,steps] = auc(a,b,plt,n)
%returns the area under the curve by sweeping a threshold through the min
%and max values of the entire dataset.  where a is the model to test and b
%is the random model to discriminate from.  A score 0f .5 is a model that
%cannot discriminate from the random distrobution.  

steps = .1;

if (nargin < 3)
    plt = 0;
end
if (nargin < 4)
    n = 'b';
end


mx = max([a;b]);
R = [];
for (ii = 0:steps:mx)   
    %tp = find((a >= ii) & (b < ii));   
    %fp = find((a >= ii) & (b >= ii));    
    %fp = find( (b >= ii) & ~((b >= ii) & (a >= ii)) );
    
    tp = find((a >= ii));
    fp = find((b >= ii));
    R = [R;[length(fp)./length(a) length(tp)./length(a)]];
end
R = [R;[0 0]];
ac = trapz(flipdim(R(:,1),1),flipdim(R(:,2),1));
steps = 0:steps:mx;


if (plt)
    plot(R(:,1),R(:,2), '-o');
    hold on;
    plot([0, 1], [0, 1], 'k--')
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
    xlim([0.0, 1.0])
    ylim([0.0, 1.0])
end
hold off;