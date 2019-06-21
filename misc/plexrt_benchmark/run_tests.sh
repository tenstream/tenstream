for dx in 100 500
do
  bash run_solver_tests.sh thermal_diff False True $dx
  bash run_solver_tests.sh solar_dir    True False $dx
  bash run_solver_tests.sh solar_diff   True False $dx
done
wait
