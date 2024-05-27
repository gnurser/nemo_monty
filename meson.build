project('nemo_monty', 'c',
  version : '0.1',
  license: 'BSD-3',
  meson_version: '>= 1.1.0',
  default_options : ['warning_level=2',
                     'buildtype=release'],
)

add_languages('fortran')
ff = meson.get_compiler('fortran')
# if ff.has_argument('-Wno-conversion')
#   add_project_arguments('-Wno-conversion', language: 'fortran')
# endif

# add_project_arguments('-march':'native', language: 'fortran')
# add_project_arguments('-fopenmp', language: 'fortran')

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()


# incdir_f2py = run_command(py,
#     ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
#     check : true
# ).stdout().strip()

#incdir_numpy = meson.get_external_property('numpy-include-dir')#, 'not-given')
incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()
incdir_f2py = incdir_numpy / '..' / '..' / 'f2py' / 'src'

inc_np = include_directories(incdir_numpy, incdir_f2py)

fortranobject_c = incdir_f2py / 'fortranobject.c'

interp_source = custom_target('interpmodule.c',
  input : ['src/nemo_monty/fortran/interp.F90'],
  output : ['interpmodule.c', 'interp-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'interp', '--lower']
)
bp_source = custom_target('bpmodule.c',
  input : ['src/nemo_monty/fortran/bp.F90'],
  output : ['bpmodule.c', 'bp-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'bp', '--lower']
)
py.extension_module('bp',
  [
    'src/nemo_monty/fortran/bp.F90', bp_source
  ],
  fortranobject_c,
  include_directories: inc_np,
  dependencies : py_dep,
  link_args : ['-lgomp'],
  install : true,
  # install_rpath: '/Users/agn/miniforge3/envs/big/lib',
  subdir: 'nemo_monty'
)

py.extension_module('interp',
  [
    'src/nemo_monty/fortran/interp.F90', interp_source
  ],
  fortranobject_c,
  include_directories: inc_np,
  dependencies : py_dep,
  link_args : ['-lgomp'],
  install : true,
  # install_rpath: '/Users/agn/miniforge3/envs/big/lib',
  subdir: 'nemo_monty'
)

python_sources = [
  'src/nemo_monty/__init__.py', 'src/nemo_monty/3dfiles_M_mp.py',
  'src/nemo_monty/findnemo.py', 'src/nemo_monty/gridstuff.py']
py.install_sources(python_sources, subdir: 'nemo_monty')

installdir = py.get_install_dir()

install_symlink('3dfiles_mp_M.py',install_dir: installdir / '../../../bin'  ,pointing_to: installdir /'nemo_monty/3dfiles_M_mp.py')

#-mtune=native -march=native
#-Wall -Wextra -Wimplicit-interface -fPIC -O3 -march=native -ffast-math -funroll-loops
