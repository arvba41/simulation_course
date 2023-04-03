clear all
% This script changes all interpreters from tex to latex. 
list_factory = fieldnames(get(groot,'factory'));

index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

f = figure(1);

plot(randn([1 100]),'DisplayName','rand')

textwidth = 12.1;
golden_ratio = (1 + sqrt(5)) / 2;
figsize = [textwidth, textwidth * golden_ratio];

% Set size and no crop
set(f, 'PaperUnits', 'centimeters', 'PaperSize', figsize);
set(f, 'PaperUnits', 'normalized', 'PaperPosition', [0, 0, 1, 1]);

print -dpdf sin_plot3.pdf