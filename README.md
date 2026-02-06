# diff-lapack

## Installation

We support building `diffblas` with the [Meson Build system](https://mesonbuild.com):

```bash
meson setup builddir
meson compile -C builddir
meson install -C builddir
```

By default, Meson compiles a static libary `libdiffblas.a`.
For a shared library, please use the option `-Ddefault_library=shared`.
To install in a specific location, use the option `--prefix=install/dir`

```bash
meson setup builddir -Ddefault_library=shared --prefix=$(pwd)/diffblas
meson compile -C builddir
meson install -C builddir
```

## Launching the differentiation

Get the latest version from the LAPACK webpage

```shell
wget -q -O tmp.html https://raw.githubusercontent.com/Reference-LAPACK/lapack/refs/heads/master/README.md
VERSION=`cat tmp.html | grep "VERSION" | tail -1 | awk '{print $3}'`
rm tmp.html
```
And use this to create the link to the most recent LAPACK release

```shell 
wget -P tmp/ https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v$VERSION.tar.gz
tar xvzf tmp/v$VERSION.tar.gz -C tmp/
rm tmp/v$VERSION.tar.gz
python run_tapenade_blas.py --input-dir tmp/lapack-$VERSION/BLAS/SRC/ --out-dir out
```

These commands are gathered in a bash file `get_and_run.sh` calling it directly from bash

Some remarks:
-------------

* If one needs to run only a subset of the BLAS files: 
```shell
python run_tapenade_blas.py --input-dir tmp/lapack-$VERSION/BLAS/SRC/ --out-dir out --file FILE1 FILE2
```
* The script has been tested for Python3. One could need to adapt the command line with the python version depending on your system.

## Calling the tests

```shell
cd out/
make # Will compile everything
./run_tests.sh # will run all available tests
```

Similarly, one can limit to certain functions only: 
```shell
make dgemm # Will compile only dgemm
```

This library generates a differentiated BLAS code using Tapenade in the following modes. 
1. Forward mode
2. Vector forward mode
3. Reverse mode
4. Vector reverse mode
