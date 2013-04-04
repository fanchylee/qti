require 'mkmf'
have_library('schmidt_decomposition', 'schmidt_decomposition')
# The destination
dir_config('qti')
# Do the work
create_makefile('qti')
