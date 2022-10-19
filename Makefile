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

# Here we create the autocorrelation figures
figures/acf_figures: .created-dirs celltrack_data/gel_data celltrack_data/glass_data acf_functions.py compile_data_tracks.py\
 libraries
	python3 acf_functions.py