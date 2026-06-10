talk: shapeconst_talk.rmd redeye_odo.Rout reedfrog.Rout
	Rscript -e "rmarkdown::render('shapeconst_talk.rmd')"

redeye_odo.Rout: odo_semimech.Rout redeye_odo.R

%.Rout: %.R
	R CMD BATCH --vanilla $<

%.html: %.rmd
	Rscript -e "rmarkdown::render('$<')"

clean:
	rm -f *~ \#*
