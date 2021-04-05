Supervised.html: Supervised.Rmd
	Rscript -e 'rmarkdown::render("Supervised.Rmd")'
Supervised-ex.html: Supervised-ex.Rmd
	Rscript -e 'rmarkdown::render("Supervised-ex.Rmd")'
depmixS4.models.RData: depmixS4.models.R
	R --no-save < $<
figure-postCP.pdf: figure-postCP.R
	R --no-save < $<
pelt.fpop.RData: pelt.fpop.R Segmentor.models.RData
	R --no-save < $<
figure-regression-interactive/index.html: figure-regression-interactive.R
	R --no-save < $<
figure-cv.png: figure-cv.R test.error.RData
	R --no-save < $<
test.error.RData: test.error.R Segmentor.models.RData
	R --no-save < $<
test.error.10fold.RData: test.error.10fold.R Segmentor.models.RData
	R --no-save < $<
Segmentor.models.RData: Segmentor.models.R # about 20 minutes on my desktop.
	R --no-save < $<
Submission.html: Submission.md
	Rscript -e 'rmarkdown::render("Submission.md")'
