function quickExport(fname)
    % Check if fname is provided
    if nargin < 1
        error('Please provide a file name as an input argument.');
    end

    % Ensure fname does not contain invalid characters for file names
    invalidChars = ['<', '>', ':', '"', '/', '\', '|', '?', '*'];
    if any(ismember(fname, invalidChars))
        error('File name contains invalid characters.');
    end

    % Set up figure for vector export
    set(gcf, 'Renderer', 'painters'); % Ensure vector rendering

    % Construct file paths dynamically
    basePath = '\\nas01.itap.purdue.edu\puhome\desktop\Figures\';
    pdfPath = fullfile(basePath, [fname, '.pdf']);
    epsPath = fullfile(basePath, [fname, '.eps']);

    % Export as PDF (Post-2022 method)
    exportgraphics(gcf, pdfPath, 'ContentType', 'vector');

    % Legacy EPS export (pre-2022)
    print('-depsc', '-painters', epsPath);

    % Notify user of successful export
    fprintf('Figure exported successfully to:\n');
    fprintf('PDF: %s\n', pdfPath);
    fprintf('EPS: %s\n', epsPath);
end
