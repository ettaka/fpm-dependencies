LIBNAME=libfpm_dependencies
PLATFORM=linux
BUILD_DIR=.
STATIC=$(LIBNAME).a ../lapack-fpm/lib/liblapack.a ../blas-fpm/lib/libblas.a ../linpack-fpm/lib/liblinpack.a

all: install shared

install: install_blas install_lapack install install_linpack install_fpm_dependencies

install_fpm_dependencies:
	fpm install --link-flag '-Llapack-fpm/lib/ -Llinpack-fpm/lib/ -Lblas-fpm/lib/' --prefix . --profile release

install_blas:
	cd blas-fpm; fpm install --prefix . --profile=release

install_lapack:
	cd lapack-fpm; fpm install --prefix . --profile=release

install_linpack:
	cd linpack-fpm; fpm install --prefix . --profile=release

shared: shared_$(PLATFORM)

shared_linux:
	cd lib; gfortran -shared -o $(LIBNAME).so -Wl,--whole-archive $(STATIC) -Wl,--no-whole-archive

shared_darwin: 
	cd lib; gfortran -dynamiclib -install_name @rpath/$(LIBNAME).dylib -static-libgfortran -static-libquadmath -static-libgcc -o $(LIBNAME).dylib -Wl,-all_load $(STATIC) -Wl,-noall_load

shared_windows: 
	cd lib; gfortran -shared -static -o $(LIBNAME).dll -Wl,--out-implib=$(LIBNAME).dll.a,--export-all-symbols,--enable-auto-import,--whole-archive $(STATIC) -Wl,--no-whole-archive

