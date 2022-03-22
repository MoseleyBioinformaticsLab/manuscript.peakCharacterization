get_assignments = function(assignments_file){
  #assignments_file = tar_read(assign_filtersd_1ecf)
  assign_data = smirfeTools::read_smirfe_assignment(assignments_file$file)
  assign_data
}
