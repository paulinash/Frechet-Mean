function colors = make_colors(N, useColor)
%MAKE_COLORS Convenience color palette for multiple curves.

if nargin < 2
    useColor = true;
end

if useColor
    % Colorful sample curves
    colors = lines(N);
    colors = 0.5*colors + 0.5; % mix with white
else
    % Gray sample curves
    gray = [0.75 0.75 0.75];
    colors = repmat(gray, N, 1);
end
end
