Segmentor.models.RData: Segmentor.models.R # about 20 minutes on my desktop.
	R --no-save < $<
Submission.html: Submission.md
	Rscript -e 'rmarkdown::render("Submission.md")'
