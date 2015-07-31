function addToPdf( dirname, pdfname )

d = dir( dirname );

input = {};
for ii = 1:length(d)
    input = {input, d(ii).name};
end

append_pdfs('out.pdf', d.name );
end