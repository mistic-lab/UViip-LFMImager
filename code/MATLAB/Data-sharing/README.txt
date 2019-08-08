## README
When condensed_data_loader.m is run in Matlab, your workspace will be
filled with ~6 variables. They are:

correlation_time    :   - the time vector from 0/fs to endtime/fs of the
                          experiment.
                        - useful for plotting against sweep_correlation

plot_magnitudes     :   - magnitudes (dB) of each correlation peak
                        - used in the colourful plot which overlays all the
                          peaks
                        - not complex

plot_xtime          :   - normalized time vector to go with plot_magnitudes
                        - each cell is 1/fs long
                        - generally from groundwave-1ms to groundwave+4ms

received_sweep      :   - raw data from the received spectrum
                        - first column is real
                        - second column is imaginary

sweep_correlation   :   - Rx and Tx files correlated
                        - complex
                        - not in dB

transmitted_sweep   :   - raw transmission for correlation
                        - first column is real
                        - second column is imaginary