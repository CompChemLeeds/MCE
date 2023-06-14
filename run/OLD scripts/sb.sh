#! /bin/bash

#Case=( 1a 1b 2 3 4 5a 5b 6a 6b 6c 7 8a 8b )
#reps=( 2000 50 50 800 50 50 50 40 40 40 50 50 50 )

#Case=( 1a 3 4 5a 5b 6a 6b 6c 7 )
#reps=( 2000 800 50 50 50 50 50 50 50 )
#spacing=( 10 30 20 30 30 20 20 10 20 )

Case=( 7 8a 8b )
reps=( 50 50 50 )
spacing=( 20 20 20 )

for i in "${!Case[@]}"; do
  cp ../SB_inputs/inputC${Case[$i]}.dat ./input.dat
  cp ../SB_inputs/inhamC${Case[$i]}.dat ./inham.dat
  cp ../SB_inputs/propC${Case[$i]}.dat ./prop.dat
  ./traintester.sh ${spacing[$i]}
  for j in 1; do
    cp ./calibinputs/input.$j input.dat
    ./run.sh ${reps[$i]} 5 5
    cp result.sh result_C${Case[$i]}_${j}.sh
  done
done
