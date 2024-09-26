#!/bin/bash

for i in `cat doc.dat`
 
do
  echo -e "\e[0;31m======================================================================================================================\e[0m"
  echo -e "\e[0;31m                                       ENTERING LAMBDA WINDOW " ${i} Please WAIT...."\e[0m"
  echo -e "\e[0;31m======================================================================================================================\e[0m"
  mkdir analysis$i
  cp ../data$i/$i.out analysis$i/ti001.out	
  cd analysis$i
  sed -n '/Summary of dvdl values over /,/End of dvdl summary/p' ti001.out >out
  sed '1d;$d' out >out1
  mv out1 ti001.out
  rm out 
  rm COLVAR
  cp ../../data$i/COLVAR .
  cp ../../data$i/HILLS HILLS
  sed -i '/^#/d' COLVAR
#  sed -i '5000001,5000050d' COLVAR                                       # Use this line when restarting the simulation.
#  sed -i '10000001,10000050d' COLVAR                                       # Use this line when restarting the simulation.
#  sed -i '15000001,15000050d' COLVAR                                       # Use this line when restarting the simulation.
  sed -i '/^#/d' HILLS 
  sed -n '1p' ti001.out > new && sed -e '1d' ti001.out | sed -n '0~10p' >>new
  mv new ti001.out 
#  cat ../../../50+50+50ns/tiufed_rew/analysis$i/ti001.out ti001.out > new      # Use these lines when restarting the simulation.
#  mv new ti001.out                                                             # Use this line when restarting the simulation.
#  sed -i '15000001,15000050d' ti001.out                                        # Use this line when restarting the simulation.
  cp -r ../3D_analysis/* . 
  cd ct_factor
  rm COLVAR HILLS
  ln -s ../COLVAR
  ln -s ../HILLS
  gfortran ct_factor.f90 -o ct_factor.x
  ./ct_factor.x
  cd ../probability/
  rm COLVAR HILLS ti001.out data_ct.dat
  ln -s ../COLVAR
  ln -s ../HILLS
  ln -s ../ti001.out
  ln -s ../ct_factor/data_ct.dat
  gfortran probability.f90 -o probability.x
  ./probability.x
  sed -i "s/XXXX/$i/g" 2_dhdl.f90 
  gfortran 2_dhdl.f90 -o 2_dhdl.x
  ./2_dhdl.x >> ../../int_dhdl.dat
  cd ../../
  echo "            lambda $i done       "
done  
gfortran deltaf.f90
./a.out
