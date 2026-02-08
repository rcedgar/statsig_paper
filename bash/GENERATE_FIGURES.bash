rm -rf ../figures
mkdir -p ../figures
python ../py/tm_gumbel_figure.py
python ../py/plot_foldseek_selected_doms.py
python ../py/3di_figure.py
python ../py/3di_frequency_stacked_hist_v2.py
python ../py/fit_logistic.py
python ../py/plot_3di_selected_scop_classes.py
python ../py/plot_tm_cmdline_doms_scop_classes.py
python ../py/plot_foldseek_cmdline_doms_scop_classes.py
python ../py/fit_tail_and_bulk_all.py
python ../py/foldseek_bits_model_ambig_tps.py
python ../py/blastp_reseek_evalue_scatterplot.py
python ../py/tm_gevd_figure_bad_scipy_fit.py
python ../py/fit_gumbel_fatcat.py
ls -lh ../figures
