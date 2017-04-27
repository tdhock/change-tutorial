Supervised.html: Supervised.Rmd
	Rscript -e 'rmarkdown::render("Supervised.Rmd")'
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
