function [ac,R] = getauc(hum, ran, name)

% loop over the random samples:
for ii = 1:size(ran,2)
  [ac(ii), R{ii}] = auc(hum, ran(:, ii));
end
