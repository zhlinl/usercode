#systematic variations (1S, unpolarized)
root -l -b -q makeWeights.C\(1,0\)     >& log/weights_AccstaLo    &
root -l -b -q makeWeights.C\(1,1\)     >& log/weights_accStaHi    &
root -l -b -q makeWeights.C\(2,0\)     >& log/weights_EtrkStaLo   &
root -l -b -q makeWeights.C\(2,1\)     >& log/weights_EtrkStaHi   &
#root -l -b -q makeWeights.C\(3,0\)     >& log/weights_EmuidStaLo  &
#root -l -b -q makeWeights.C\(3,1\)     >& log/weights_EmuidStaHi  &
#root -l -b -q makeWeights.C\(4,0\)     >& log/weights_EtrigstaLo  &
#root -l -b -q makeWeights.C\(4,1\)     >& log/weights_EtrigstaHi  &
root -l -b -q makeWeights.C\(5,0\)     >& log/weights_EtrecostaLo  &
root -l -b -q makeWeights.C\(5,1\)     >& log/weights_EtrecostaHi  &
root -l -b -q makeWeights.C\(6,0\)     >& log/weights_PtscaleLo   &
root -l -b -q makeWeights.C\(6,1\)     >& log/weights_PtscaleHi   &
root -l -b -q makeWeights.C\(7,0\)     >& log/weights_PtresoLo    &
root -l -b -q makeWeights.C\(7,1\)     >& log/weights_PtresoHi    &
root -l -b -q makeWeights.C\(8,0\)     >& log/weights_Ptspec      &
root -l -b -q makeWeights.C\(9,0\)     >& log/weights_VtxPos      &
root -l -b -q makeWeights.C\(10,0\)    >& log/weights_Nofsr	      &
root -l -b -q makeWeights.C\(11,0\)    >& log/weights_Tnpmc	      &
root -l -b -q makeWeights.C\(12,0\)    >& log/weights_Mctrue      &
root -l -b -q makeWeights.C\(13,0\)    >& log/weights_Tnpmcups    &
root -l -b -q makeWeights.C\(14,0\)    >& log/weights_bgLinear    &
root -l -b -q makeWeights.C\(15,0\)    >& log/weights_otherLo     &
root -l -b -q makeWeights.C\(15,1\)    >& log/weights_otherHi     &
#1s					
root -l -b -q makeWeights.C\(0,0,0,0\) >& log/weights_unpol_1s    &
root -l -b -q makeWeights.C\(0,0,0,1\) >& log/weights_helT_1s     &
root -l -b -q makeWeights.C\(0,0,0,2\) >& log/weights_helL_1s     &
root -l -b -q makeWeights.C\(0,0,0,3\) >& log/weights_csT_1s      &
root -l -b -q makeWeights.C\(0,0,0,4\) >& log/weights_csL_1s      &
#2s							      
root -l -b -q makeWeights.C\(0,0,1,0\) >& log/weights_unpol_2s    &
root -l -b -q makeWeights.C\(0,0,1,1\) >& log/weights_helT_2s     &
root -l -b -q makeWeights.C\(0,0,1,2\) >& log/weights_helL_2s     &
root -l -b -q makeWeights.C\(0,0,1,3\) >& log/weights_csT_2s      &
root -l -b -q makeWeights.C\(0,0,1,4\) >& log/weights_csL_2s      &
#3s							      
root -l -b -q makeWeights.C\(0,0,2,0\) >& log/weights_unpol_3s    &
root -l -b -q makeWeights.C\(0,0,2,1\) >& log/weights_helT_3s     &
root -l -b -q makeWeights.C\(0,0,2,2\) >& log/weights_helL_3s     &
root -l -b -q makeWeights.C\(0,0,2,3\) >& log/weights_csT_3s      &
root -l -b -q makeWeights.C\(0,0,2,4\) >& log/weights_csL_3s      &  

