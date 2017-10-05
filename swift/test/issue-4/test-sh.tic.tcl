if { $env(PMI_RANK) == 0 } {
  exec ./g09.sh mpirun.inp
}
