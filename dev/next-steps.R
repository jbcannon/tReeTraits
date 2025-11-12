## NExt steps for QSM
# - address singular matrix issue w/ hagood
# - build help for multiple optimization types
# output/assess fit for multiple parameters

install_PyTLidar()
#install_PyTLidar(reinstall=TRUE)

file = 'inst/extdata/tree_0129.laz'
ipd = NULL; pdmin = NULL; pdmax = NULL; brad1 = NULL; brad2 = NULL; plot = NULL

x = run_treeqsm(file, patch_diam2min=0.04, patch_diam2max = 0.12)

# This writes the QSM reconstruction, parameter table to disk.
write_qsm(x, name = 'testTreeID',output_dir = output_dir)

treeID = 'testTreeID'


write_qsm(x$qsm, name='test', output_dir)


load_qsm = function(path) {
  qsm = readr::read_delim(path)

}
