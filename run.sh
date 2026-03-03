#! /bin/bash

vpkg_require gcc/4.9
#vpkg_require python/20221031:openqs
#gfortran -o run_phi_program main_program.f90
#echo lokooz2

vpkg_devrequire intel
export FC=ifort
export FFLAGS=-mkl


ipara=$(awk '/iparallel/{print $1}' spectra_cmrt.inp)
iflow=$(awk '/iflow/{print $1}' spectra_cmrt.inp)
imaxproc=$(awk '/iproc/{print $1}' spectra_cmrt.inp)
iwait=$(awk '/iwait/{print $1}' spectra_cmrt.inp)

crdir()
{
  for((i=1;i<=$ipara;i++)); do
    if [ -d $i ]; then
      rm -r $i
      mkdir $i
    else
      mkdir $i
    fi
    cd $i
    cp ../aout .
    cp ../*.inp .
    cp ../submit.qs .
    echo $i > ifolder.inp
    sbatch submit.qs
    sleep 0.5
    cd ..
  done
}


avg()
{
  if [ $ipara == 1 ]; then
    echo "Cannot average. Not parallel job",$ipara
  else
    rm response.out se.out gsb.out esa.out
    for((i=1;i<=$ipara;i++)); do
      cat $i/response.out >> response.out
      cat $i/se.out >> se.out
      cat $i/gsb.out >> gsb.out
      cat $i/esa.out >> esa.out
      #file_ra=$file_ra" $i/pop.out"
    done
  fi

}

variation_param()
{
  awk 'NR==1{$1=2}1' spectra_cmrt.inp>save_spectra_cmrt.inp

  while read lineno newvar dire
    do
      mkdir -p $dire
      cd $dire
      cp ../run .
      cp ../submit.qs .
      cp ../*.f* .
      cp ../makefile .
      touch mod*
      awk 'NR==n{$1=x}1' n=$lineno x=$newvar ../save_spectra_cmrt.inp>spectra_cmrt.inp
      ./run
      cd ..
    done < variation.inp
}

if [ $iflow == 1  ]; then
  ### serial execution
  echo "serial execution"
  make
  sbatch submit.qs
elif [ $iflow == 2 ];then
  ### parallel execution
  echo "parallel execution, Number of cpus = "$ipara
  make
  crdir
elif [ $iflow == 3 ]; then
  ### average for parallel execution
  echo "Averaging data"
  avg
elif [ $iflow == 4 ]; then
  ### variation to parameters
  echo "Variation to parameters"
  variation_param
else
  echo "no matching iflow",$iflow
fi

