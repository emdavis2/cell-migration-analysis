.PHONY: clean
.PHONY: d3-vis
.PHONY: visualization

clean:
	rm -rf model
	rm -rf hetero_model
	rm -rf figures
	rm -rf .created-dirs
	rm -f writeup.aux
	rm -f writeup.log
	rm -f writeup.pdf

.created-dirs:
	mkdir -p model
	mkdir -p hetero_model
	mkdir -p figures
	mkdir -p figures/acf_figures
	mkdir -p figures/histogram_boxplot
	mkdir -p figures/model
	mkdir -p figures/hetero_model
	touch .created-dirs

# Create the autocorrelation figures for glass data
figures/acf_figures_glass: .created-dirs celltrack_data/glass_data\
 functions/acf_functions.py functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py
	python3.7 GenerateDataACF.py 'celltrack_data/glass_data' 30 'glass'

# Create the autocorrelation figures for gel data
figures/acf_figures_gel: .created-dirs celltrack_data/gel_data\
 functions/acf_functions.py functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py
	python3.7 GenerateDataACF.py 'celltrack_data/gel_data' 30 'stiff'

# Create the boxplot and histogram figures for both glass and gel data
figures/histogram_boxplot: .created-dirs celltrack_data/gel_data\
 celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py
	python3.7 GenerateDataHistogramBoxplot.py 'celltrack_data/glass_data' 'celltrack_data/gel_data' 30 'glass' 'stiff'

#Perform grid search to fit PRW model to glass data 
model/vel_acf_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW' 'run_sim_get_velacf_err'

#Perform grid search to fit PRW model to gel data 
model/vel_acf_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW' 'run_sim_get_velacf_err'

#Perform grid search to fit PRW_polaritybias model to glass data 
model/vel_acf_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'run_sim_get_velacf_err'

#Perform grid search to fit PRW_polaritybias model to gel data 
model/vel_acf_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'run_sim_get_velacf_err'

#Make figures compring models with optimal parameters to glass data
figures/model_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_PRW.txt model/model_params_glass_PRW_PB.txt
	python3.7 CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'model/model_params_glass_PRW.txt' 'model/model_params_glass_PRW_PB.txt'

#Make figures compring models with optimal parameters to gel data
figures/model_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_PRW.txt model/model_params_stiff_PRW_PB.txt
	python3.7 CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'model/model_params_stiff_PRW.txt' 'model/model_params_stiff_PRW_PB.txt'

# Build the final report for the project.
writeup.pdf: figures/acf_figures/glass_polarity_vector_acf_avg.png figures/acf_figures/stiff_polarity_vector_acf_avg.png\
 figures/acf_figures/glass_polarity_angle_acf_avg.png figures/acf_figures/stiff_polarity_angle_acf_avg.png\
 figures/acf_figures/glass_abs_skew_acf_avg.png figures/acf_figures/stiff_abs_skew_acf_avg.png\
 figures/acf_figures/glass_velocity_angle_acf_avg.png figures/acf_figures/stiff_velocity_angle_acf_avg.png\
 figures/acf_figures/glass_velocity_acf_avg.png figures/acf_figures/stiff_velocity_acf_avg.png\
 figures/acf_figures/glass_speed_acf_avg.png figures/acf_figures/stiff_speed_acf_avg.png\
 figures/acf_figures/glass_speed_x_acf_avg.png figures/acf_figures/stiff_speed_acf_avg.png\
 figures/acf_figures/glass_speed_y_acf_avg.png figures/acf_figures/stiff_speed_y_acf_avg.png
	pdflatex writeup.tex

#Perform grid search to fit heterogeneous PRW model to glass data 
hetero_model/vel_acf_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW' 'run_sim_get_velacf_err'

#Perform grid search to fit heterogeneous PRW model to gel data 
hetero_model/vel_acf_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW' 'run_sim_get_velacf_err'

#Perform grid search to fit heterogeneous PRW_polaritybias model to glass data 
hetero_model/vel_acf_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'run_sim_get_velacf_err'

#Perform grid search to fit heterogeneous PRW_polaritybias model to gel data 
hetero_model/vel_acf_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'run_sim_get_velacf_err'

#Make figures compring heterogeneous models with optimal parameters to glass data
figures/hetero_model_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_PRW.txt model/model_params_glass_PRW_PB.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'hetero_model/model_params_glass_PRW.txt' 'hetero_model/model_params_glass_PRW_PB.txt'

#Make figures compring heterogeneous models with optimal parameters to gel data
figures/hetero_model_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_PRW.txt model/model_params_stiff_PRW_PB.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'hetero_model/model_params_stiff_PRW.txt' 'hetero_model/model_params_stiff_PRW_PB.txt'