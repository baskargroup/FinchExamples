function cb = add_colorbar(ax, S, label)
% Colour bar in the shared style.
    cb = colorbar(ax); cb.Label.String = label;
    cb.Label.Color = S.ink; cb.Color = S.ink2; cb.TickDirection = 'in';
    cb.FontSize = S.fs;
end
