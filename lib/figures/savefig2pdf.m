function savefig2pdf(fig_handle, fname)
% save figure to pdf with appropriate resizing of paper to eliminate whitespace    

if ~isvalid(fig_handle)
    error('fig_handle to savefig2pdf.m was not valid');
end

% TODO: error if not .pdf or empty filename

fig_handle.Units = 'Inches';
fig_handle.PaperPositionMode = 'Auto';
fig_handle.PaperUnits = 'Inches';
fig_handle.PaperSize = fig_handle.Position(3:4);

print(fig_handle, fname, '-dpdf', '-r0');

end