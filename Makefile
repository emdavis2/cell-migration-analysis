.PHONY: clean
.PHONY: d3-vis
.PHONY: visualization

clean:
	rm -rf model
	rm -rf figures
	rm -rf derived_data
	rm -rf .created-dirs
	rm -f writeup.pdf

.created-dirs:
	mkdir -p model
	mkdir -p figures
	mkdir -p derived_data
	touch .created-dirs

# It is usefual to propogate the panas_na and panas_pa values forward
# in time for subsequent analysis which attempts to cluster the state
# of each person in the study regardless of time so that we can track
# their progress through a lower dimensional space.
derived_data/clinical-outcomes-preprocessed.csv: .created-dirs pre-process.R source_data/clinical_outcomes.csv
	Rscript pre-process.R