function lettering(letter)
x_shift = -0.08;
y_shift = 0.11;


annotation('textbox', get(gca, 'Position') + [x_shift y_shift 0 0],...
    'String', letter, 'FitBoxToText', 'on',...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top',...
    'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold',...
    'LineStyle', 'none');