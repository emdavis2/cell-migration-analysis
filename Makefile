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
	mkdir -p figures/velacf_model
	mkdir -p figures/MSD_model
	mkdir -p figures/velacf_hetero_model
	mkdir -p figures/MSD_hetero_model
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

#Perform grid search to fit PRW model to glass data with vel acf
model/vel_acf_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'LPRW' 'vel_acf' 'S' 10 50 20 'P' 0.5 5 20

#Perform grid search to fit PRW model to gel data with vel acf
model/vel_acf_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'LPRW' 'vel_acf' 'S' 10 50 20 'P' 0.5 5 20

#Perform grid search to fit PRW_polaritybias model to glass data with vel acf
model/vel_acf_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'vel_acf' 'std_dev_w' 0.2 0.9 10 'std_dev_theta' 0.9 1.5 10

#Perform grid search to fit PRW_polaritybias model to gel data with vel acf
model/vel_acf_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'vel_acf' 'std_dev_w' 0.2 0.9 10 'std_dev_theta' 0.9 1.5 10

#Perform grid search to fit PRW model to glass data with MSD
model/MSD_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'LPRW' 'MSD' 'S' 10 50 20 'P' 0.5 5 20

#Perform grid search to fit PRW model to gel data with MSD
model/MSD_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'LPRW' 'MSD' 'S' 10 50 20 'P' 0.5 5 20

#Perform grid search to fit PRW_polaritybias model to glass data with MSD
model/MSD_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'MSD' 'std_dev_w' 0.2 0.9 10 'std_dev_theta' 0.9 1.5 10

#Perform grid search to fit PRW_polaritybias model to gel data with MSD
model/MSD_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'MSD' 'std_dev_w' 0.2 0.9 10 'std_dev_theta' 0.9 1.5 10

#Make figures compring models with optimal parameters to glass data using vel_acf fitting metric
figures/model_glass_velacf: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_LPRW_vel_acf.txt model/model_params_glass_PRW_PB_vel_acf.txt
	python3.7 CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'model/model_params_glass_LPRW_vel_acf.txt' 'model/model_params_glass_PRW_PB_vel_acf.txt' 'figures/velacf_model/'

#Make figures compring models with optimal parameters to gel data using vel_acf fitting metric
figures/model_gel_velacf: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_LPRW_vel_acf.txt model/model_params_stiff_PRW_PB_vel_acf.txt
	python3.7 CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'model/model_params_stiff_LPRW_vel_acf.txt' 'model/model_params_stiff_PRW_PB_vel_acf.txt' 'figures/velacf_model/'

#Make figures compring models with optimal parameters to glass data using MSD fitting metric
figures/model_glass_MSD: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_LPRW_MSD.txt model/model_params_glass_PRW_PB_MSD.txt
	python3.7 CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'model/model_params_glass_LPRW_MSD.txt' 'model/model_params_glass_PRW_PB_MSD.txt' 'figures/MSD_model/'

#Make figures compring models with optimal parameters to gel data using MSD fitting metric
figures/model_gel_MSD: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_LPRW_MSD.txt model/model_params_stiff_PRW_PB_MSD.txt
	python3.7 CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'model/model_params_stiff_LPRW_MSD.txt' 'model/model_params_stiff_PRW_PB_MSD.txt' 'figures/MSD_model/'

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

#Perform grid search to fit heterogeneous PRW model to glass data with vel acf
hetero_model/vel_acf_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'LPRW' 'vel_acf' 'S' 10 50 10 'P' 0.5 5 10

#Perform grid search to fit heterogeneous PRW model to gel data with vel acf
hetero_model/vel_acf_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'LPRW' 'vel_acf' 'S' 10 50 10 'P' 0.5 5 10

#Perform grid search to fit heterogeneous PRW_polaritybias model to glass data with vel acf
hetero_model/vel_acf_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'vel_acf' 'std_dev_w' 0.1 2 10 'std_dev_theta' 0.1 2 10

#Perform grid search to fit heterogeneous PRW_polaritybias model to gel data with vel acf
hetero_model/vel_acf_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'vel_acf' 'std_dev_w' 0.1 2 10 'std_dev_theta' 0.1 2 10

#Perform grid search to fit heterogeneous PRW model to glass data with MSD
hetero_model/MSD_PRW_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW' 'MSD' 'S' 10 50 10 'P' 0.5 5 10

#Perform grid search to fit heterogeneous PRW model to gel data with MSD
hetero_model/MSD_PRW_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/hetero_model_fitting_functions.py functions/PRW_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW' 'MSD' 'S' 10 50 10 'P' 0.5 5 10

#Perform grid search to fit heterogeneous PRW_polaritybias model to glass data with MSD
hetero_model/MSD_PRW_PB_glass: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'PRW_PB' 'MSD' 'std_dev_w' 0.1 2 10 'std_dev_theta' 0.1 2 10

#Perform grid search to fit heterogeneous PRW_polaritybias model to gel data with MSD
hetero_model/MSD_PRW_PB_gel: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/msd_functions.py functions/hetero_model_fitting_functions.py functions/PRWpolaritybias_model_functions.py
	python3.7 Hetero_RunGridSearchFitModel.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'PRW_PB' 'MSD' 'std_dev_w' 0.1 2 10 'std_dev_theta' 0.1 2 10

#Make figures compring heterogeneous models with optimal parameters to glass data using vel_acf fitting metric
figures/hetero_model_glass_velacf: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_PRW_vel_acf.txt model/model_params_glass_PRW_PB_vel_acf.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'hetero_model/model_params_glass_PRW_vel_acf.txt' 'hetero_model/model_params_glass_PRW_PB_vel_acf.txt' 'figures/velacf_hetero_model/'

#Make figures compring heterogeneous models with optimal parameters to gel data using vel_acf fitting metric
figures/hetero_model_gel_velacf: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_PRW_vel_acf.txt model/model_params_stiff_PRW_PB_vel_acf.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'hetero_model/model_params_stiff_PRW_vel_acf.txt' 'hetero_model/model_params_stiff_PRW_PB_vel_acf.txt' 'figures/velacf_hetero_model/'

#Make figures compring heterogeneous models with optimal parameters to glass data using MSD fitting metric
figures/hetero_model_glass_MSD: .created-dirs celltrack_data/glass_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_glass_PRW_MSD.txt model/model_params_glass_PRW_PB_MSD.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/glass_data' 30 'glass' 5 0.1667 113 'hetero_model/model_params_glass_PRW_MSD.txt' 'hetero_model/model_params_glass_PRW_PB_MSD.txt' 'figures/MSD_hetero_model/'

#Make figures compring heterogeneous models with optimal parameters to gel data using MSD fitting metric
figures/hetero_model_gel_MSD: .created-dirs celltrack_data/gel_data functions/compile_data_tracks_function.py\
 functions/libraries/track_functions.py functions/libraries/qc_functions.py\
 functions/libraries/filter_cells_fns.py functions/libraries/centers.py\
 functions/acf_functions.py functions/msd_functions.py functions/PRWpolaritybias_model_functions.py functions/PRW_model_functions.py\
 functions/model_fitting_functions.py model/model_params_stiff_PRW_MSD.txt model/model_params_stiff_PRW_PB_MSD.txt
	python3.7 Hetero_CompareModelAndData.py 'celltrack_data/gel_data' 30 'stiff' 5 0.1667 119 'hetero_model/model_params_stiff_PRW_MSD.txt' 'hetero_model/model_params_stiff_PRW_PB_MSD.txt' 'figures/MSD_hetero_model/'