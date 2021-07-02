plot_hist -i state_A/configuration.dat state_B/configuration.dat -x "Bond length (nm)" -t "Bond length 1-2" -n BL_12_hist -l "State A" "State B" -c 1 -r 0.12 0.19 -k
plot_hist -i state_A/configuration.dat state_B/configuration.dat -x "Bond length (nm)" -t "Bond length 2-3" -n BL_23_hist -l "State A" "State B" -c 2 -r 0.12 0.19 -k
plot_hist -i state_A/configuration.dat state_B/configuration.dat -x "Bond length (nm)" -t "Bond length 3-4" -n BL_34_hist -l "State A" "State B" -c 3 -r 0.11 0.19 -k
plot_hist -i state_A/configuration.dat state_B/configuration.dat -x "Angle (deg)" -cc "radian to degree" -t "Angle 1-2-3" -n angle_123_hist -l "State A" "State B" -c 4 -r 80 160 -k
plot_hist -i state_A/configuration.dat state_B/configuration.dat -x "Angle (deg)" -cc "radian to degree" -t "Angle 2-3-4" -n angle_234_hist -l "State A" "State B" -c 5 -r 80 160 -k 
combine_plots -f angle_* -b -n angle_hist
combine_plots -f BL_* -b -n bond_length_hist -d 3 1


