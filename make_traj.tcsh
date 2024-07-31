#!/usr/bin/tcsh

set output_files = `ls -1 output_*.pdb | sort -V`

set trajectory_file = "traj.pdb"
if (-e $trajectory_file) then
  rm $trajectory_file
endif

foreach file ( $output_files )
  echo "Processing $file ..."
  cat $file >> $trajectory_file
end

echo "Trajectory saved to $trajectory_file."

# rm output_*.pdb

