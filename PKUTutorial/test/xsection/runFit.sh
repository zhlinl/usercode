#pt-diff, 1 y-bin
./runFitSyst.py 1 0  >& log/fitsys_1s_mode0 &
./runFitSyst.py 2 0  >& log/fitsys_2s_mode0 &
./runFitSyst.py 3 0  >& log/fitsys_3s_mode0 &

#pt-diff, 2 y-bin
./runFitSyst.py 1 1  >& log/fitsys_1s_mode1 &
./runFitSyst.py 2 1  >& log/fitsys_2s_mode1 &
./runFitSyst.py 3 1  >& log/fitsys_3s_mode1 &

#y-diff, 1 pt bin
./runFitSyst.py 1 2  >& log/fitsys_1s_mode2 &
./runFitSyst.py 2 2  >& log/fitsys_2s_mode2 &
./runFitSyst.py 3 2  >& log/fitsys_3s_mode2 &

#ratios
./runFitSyst.py 1 3  >& log/fitsys_1s_mode3 &
./runFitSyst.py 2 3  >& log/fitsys_2s_mode3 &
./runFitSyst.py 3 3  >& log/fitsys_3s_mode3 &

#ps | grep oniafitter.C | grep root | grep -v exe|wc
#ls -1 mode1/fitres_?s_*/xsecdiff*rap0.pdf|wc #74
#ls -1 mode?/fitres_?s_*/xsecdiff*rap0.pdf|wc #296
