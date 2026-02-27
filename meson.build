project(
  'diffblas',
  'fortran', 'c',
  version: '2026.2.27',
  meson_version: '>= 0.63.0',
  default_options: [
    'buildtype=release',
    'libdir=lib',
    'default_library=static',
    'warning_level=0',
  ],
)

# Compilers
fc = meson.get_compiler('fortran')

# Recognise old non-standard double complex intrinsics
if fc.get_id() == 'nagfor'
  add_global_arguments('-dcfuns', language : 'fortran')
endif

# Options
install_modules = get_option('modules')
libblas_name = get_option('libblas')
libblas_path = get_option('libblas_path')
libcblas_name = get_option('libcblas')
libcblas_path = get_option('libcblas_path')

# Dependencies
if libblas_path == []
    libblas = dependency(libblas_name, required : false)
endif
if libblas_path != [] or not libblas.found()
    libblas = fc.find_library(libblas_name, dirs : libblas_path, required : false)
endif
if libcblas_path == []
    libcblas = dependency(libcblas_name, required : false)
endif
if libcblas_path != [] or not libcblas.found()
    libcblas = fc.find_library(libcblas_name, dirs : libcblas_path, required : false)
endif

# Headers
libdiffblas_include = include_directories('BLAS/include', 'TAPENADE/include')
libdiffcblas_include = include_directories('CBLAS/src/tmp', 'TAPENADE/include')

# Sources
libdiffblas_src = []
libdiffcblas_src = []

subdir('BLAS')
subdir('CBLAS')
subdir('TAPENADE')

# Library
libdiffblas = library('diffblas',
                      sources : libdiffblas_src,
                      dependencies : libblas,
                      include_directories: libdiffblas_include,
                      install : true)

libdiffcblas = library('diffcblas',
                      sources : libdiffcblas_src,
                      dependencies : libblas,
                      include_directories: libdiffcblas_include,
                      install : true)

# Fortran modules
if install_modules
  script_modules = files('install_modules.py')
  meson.add_install_script(script_modules)
endif
