#!/bin/bash

for j in `seq 1000000 500000 5000000`
do
  for i in `cat doc.dat`
   
  do
    echo -e "\e[0;31m======================================================================================================================\e[0m"
    echo -e "\e[0;31m                                       ENTERING LAMBDA WINDOW " ${i} Please WAIT...."\e[0m"
    echo -e "\e[0;31m======================================================================================================================\e[0m"
#    mkdir analysis$i
  #  sed '1d;$d' ../data$i/ti001.out > analysis$i/ti001.out
#    cp ../data$i/$i.out analysis$i/ti001.out	
    cd analysis$i
#    sed -n '/Summary of dvdl values over /,/End of dvdl summary/p' ti001.out >out
#    sed '1d;$d' out >out1
#    mv out1 ti001.out
#    rm out 
#    rm COLVAR
#    cp ../../data$i/COLVAR* .
#    cp ../../data$i/HILLS HILLS
#    sed -i '/^#/d' COLVAR*
#    sed -i '/^#/d' HILLS 
  #  ln -s ../../data$i/COLVAR .
#    sed -n '1p' ti001.out > new && sed -e '1d' ti001.out | sed -n '0~10p' >>new
#    mv new ti001.out 
    cp -r ../3D_analysis/probability/input_conv probability/
    cp -r ../3D_analysis/probability/2_dhdl.f90 probability/
#    cd ct_factor
#    rm COLVAR HILLS
#    ln -s ../COLVAR_CC
#    ln -s ../HILLS
#    gfortran ct_factor.f90 -o ct_factor.x
#    ./ct_factor.x
    cd probability/
    cp input_conv input
    sed -i -e "s/XXXX/$j/g" input
#    rm COLVAR HILLS ti001.out data_ct.dat
#    ln -s ../COLVAR
#    ln -s ../HILLS
#    ln -s ../ti001.out
#    ln -s ../ct_factor/data_ct.dat
    gfortran probability.f90 -o probability.x
    ./probability.x
    sed -i "s/XXXX/$i/g" 2_dhdl.f90 
    gfortran 2_dhdl.f90 -o 2_dhdl.x
    ./2_dhdl.x >> ../../int_dhdl.dat
    cd ../../
    echo "            lambda $i done       "
  done  
  cp deltaf_conv.f90 deltaf.f90
  sed -i -e "s/XXXX/$j/g" deltaf.f90
  gfortran deltaf.f90
  ./a.out >>out
   mv int_dhdl.dat int_dhdl$j.dat
   mv deltaf.dat deltaf$j.dat
done
