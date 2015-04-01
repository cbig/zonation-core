# Depends on: boost (program_options, filesystem, system, chrono (since boost 1.51), thread), QtCore, QtNetwork, gdal, fftw3/mkl

include_directories(${PROJECT_SOURCE_DIR})

file(GLOB ZIG4_SRCS *.cpp core/matrix_utils.c core/randz.c core/*.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/core)

add_executable(zig4 ${ZIG4_SRCS})

if (CMAKE_COMPILER_IS_GNUCC)
  if(NOT COMPILE_STATIC)  
    #set(CXX_OPT_OPT "-g  -march=core2 -mfpmath=sse -mieee-fp -funroll-loops -pipe -mmmx -msse -msse2 -msse3 -m3dnow")  # More risk: -mlwp -mabm -mbmi -mtbm    / i3/i5/i7: avx!
    set(CXX_OPT_OPT "-O2 -finline-functions -march=core2 -mfpmath=sse -mieee-fp -funroll-loops -pipe -mmmx -msse -msse2 -msse3 -m3dnow")  # More risk: -mlwp -mabm -mbmi -mtbm    / i3/i5/i7: avx!
  else()
    # w/o -mieee-fp
    set(CXX_OPT_OPT "-O2 -finline-functions -march=core2 -mfpmath=sse -funroll-loops -pipe -mmmx -msse -msse2 -msse3 -m3dnow")  # More risk: -mlwp -mabm -mbmi -mtbm    / i3/i5/i7: avx!
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_OPT_OPT}")  
endif()

find_package(Qt4 REQUIRED QtNetwork QtCore)
include(${QT_USE_FILE})

if (COMPACT_VMAT_LOOKUP)
  add_definitions(-DCOMPACT_VMAT)
  add_definitions(-DCOMPACT_VMAT_LOOKUP)
endif()
if (NONCOMPACT_VMAT_FLOATP)
  add_definitions(-DCOMPACT_VMAT)
  add_definitions(-DNONCOMPACT_VMAT_FLOATP)
endif()

if (NOT MINGW)
   find_package(GDAL REQUIRED)
   if(GDAL_FOUND)
     include_directories(${GDAL_INCLUDE_DIR})
   endif()
endif()

if (MINGW)
   set(Boost_ADDITIONAL_VERSIONS 1_51 1_51_0)
   set(Boost_DIR /mingw)
endif()

# Link
set(Boost_USE_MULTITHREADED ON)
if (MINGW)
   set(Boost_USE_STATIC_LIBS ON)
   set(Boost_USE_STATIC_RUNTIME ON)
else()
   set(Boost_USE_STATIC_LIBS OFF)
   set(Boost_USE_STATIC_RUNTIME OFF)
endif()
find_package(Boost 1.45.0 COMPONENTS program_options filesystem system chrono thread)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
elseif (MSYS)
  include_directories(/mingw/include)
  set(Boost_LIBRARIES ${LOCAL_LIBDIR}/libboost_program_options-mgw48-mt-1_55.a ${LOCAL_LIBDIR}/libboost_thread-mgw48-mt-1_55.a ${LOCAL_LIBDIR}/libboost_filesystem-mgw48-mt-1_55.a ${LOCAL_LIBDIR}/libboost_system-mgw48-mt-1_55.a ${LOCAL_LIBDIR}/libboost_chrono-mgw48-mt-1_55.a)
endif()

option(USE_OPENMP
	"Use OpenMP" ON)
if(USE_OPENMP)
	include(FindOpenMP)
	if(OPENMP_FOUND)
		set({CMAKE_C_FLAGS} "${CMAKE_C_FLAGS} ${OepnMP_C_FLAGS}")
		set({CMAKE_CXX_FLAGS} "${CMAKE_CXX_FLAGS} ${OepnMP_CXX_FLAGS}")
		set({CMAKE_EXE_LINKER_FLAGS} "${CMAKE_EXE_LINKER_FLAGS} ${OepnMP_EXE_LINKER_FLAGS}")
	endif()
endif()

# To build redistributable closed source binaries
# Requires a licensed MKL - which doesn't work well anyway
if (MINGW)
   set(USE_MKL_FFT 0)
endif()

if (MINGW)
   set(QT_LIBRARIES ${LOCAL_LIBDIR}/libQtNetwork.a ${LOCAL_LIBDIR}/libQtCore.a ws2_32)
   if (MINGW_IS_W64)
      if (USE_MKL_FFT)
	    # from "Intel compiler"/compiler/lib/intel64
	    # Static
	    set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl/intel64/fftw3xf_intel.lib -Wl,--start-group ${LOCAL_LIBDIR}/mkl/intel64/mkl_sequential.lib  ${LOCAL_LIBDIR}/mkl/intel64/mkl_intel_ilp64.lib ${LOCAL_LIBDIR}/mkl/intel64/mkl_core.lib ${LOCAL_LIBDIR}/mkl/intel64/chkstk.obj  -Wl,--end-group)
	    # Dynamic - dll
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_intel_ilp64_dll.lib ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_sequential_dll.lib ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_core_dll.lib) 
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl/intel64/fftw3xf_intel.lib ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_intel_ilp64_dll.lib ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_sequential_dll.lib ${LOCAL_LIBDIR}/mkl-dll/intel64/mkl_core_dll.lib) # ${LOCAL_LIBDIR}/mkl/intel64/dll/mkl_cdft_core.dll ${LOCAL_LIBDIR}/mkl/intel64/dll/mkl_cdft_core.dll) 

      else()
            set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/libfftw3f.a)
      endif()

      set(GDAL_LIBRARIES ${LOCAL_LIBDIR}/libgdal.a wsock32 ws2_32)
      target_link_libraries(zig4 zig4lib ${FFTW3_LIBRARIES} ${Boost_LIBRARIES} ${QT_LIBRARIES} ${GDAL_LIBRARIES})

   else()  # mingw 32 bits
      if (USE_MKL_FFT)
	    # compiled from "Intel compiler"/mkl/interfaces/fftw3xf/
	    # Follow the "Intel mkl library link line advisor" - online - it's imcomplete and outdated!

	    # Static (started) based on dynamic dll dependencies
	    # Sequential!   mkl_sequential.lib > mkl_intel_thread.lib
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl/ia32/fftw3xf_intel.lib -Wl,--start-group ${LOCAL_LIBDIR}/mkl/ia32/mkl_sequential.lib  ${LOCAL_LIBDIR}/mkl/ia32/mkl_intel_c.lib ${LOCAL_LIBDIR}/mkl/ia32/mkl_core.lib ${LOCAL_LIBDIR}/mkl/ia32/chkstk.obj ${LOCAL_LIBDIR}/mkl/ia32/lldiv.obj ${LOCAL_LIBDIR}/mkl/ia32/llrem.obj ${LOCAL_LIBDIR}/mkl/ia32/llmul.obj -Wl,--end-group) # ${LOCAL_LIBDIR}/mkl/ia32/msvcrt.lib )   #${LOCAL_LIBDIR}/mkl/ia32/mkl_rt.lib #mkl_intel_c mkl_sequential
	    # ***** does not work: dynamic (because of mkl_rt and mkl_intel_thread!)
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl/ia32/fftw3xf_intel.lib ${LOCAL_LIBDIR}/mkl/ia32/mkl_rt.dll ${LOCAL_LIBDIR}/mkl/ia32/mkl_intel_thread.dll)   #mkl_intel_c mkl_sequential 

	    # Dynamic - dll
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl-dll/ia32/mkl_intel_ilp64_dll.lib ${LOCAL_LIBDIR}/mkl-dll/ia32/mkl_sequential_dll.lib ${LOCAL_LIBDIR}/mkl-dll/ia32/mkl_core_dll.lib) 
	    #set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/mkl/ia32/dll/mkl_cdft_core.dll ${LOCAL_LIBDIR}/mkl/ia32/dll/mkl_rt.dll)

      else()
            set(FFTW3_LIBRARIES ${LOCAL_LIBDIR}/libfftw3f.a) # ${LOCAL_LIBDIR}/libpthreadGC2.a # Would use pthread!: fftw3f_omp
      endif()
      set(GDAL_LIBRARIES ${LOCAL_LIBDIR}/libgdal.a)
      target_link_libraries(zig4 ${FFTW3_LIBRARIES} ${Boost_LIBRARIES} ${QT_LIBRARIES} ${GDAL_LIBRARIES} ${LOCAL_LIBDIR}/libiconv.a zig4lib ${Boost_LIBRARIES}  )  # ${GDAL_LIBRARIES}
   endif()
else()
  if (COMPILE_STATIC)
    # For gdal related libs, ask: 'gdal-config --dep-libs' OR compile a tailored libgdal.a
    # example configuration for gdal: ./configure --without-sqlite3 --without-ogdi --without-java --without-armadillo  --without-xerces --without-mysql  --without-odbc --without-pam --without-opencl --without-curl  --without-spatialite --without-pg  --without-vfk --without-freexl  --without-netcdf
    #  --with-static-proj4 ?
    # ogdi turns out to be hard to compile static. To begin with, needs:
    # CFLAGS="-fPIC" ./configure --with-proj && make
    # target_link_libraries(zig4 zig4lib QtNetwork QtCore z boost_program_options boost_thread-mt boost_filesystem-mt boost_system-mt fftw3f fftw3f_threads pthread gdal -lgeos_c -lexpat -ljasper -lhdf5 -lmfhdfalt -ldfalt -lgif -ljpeg -lpng -lz -lpthread -lm -lrt -ldl -static-libgcc -static-libstdc++ )

    set(LOCAL_STATIC_LIBDIR "/home/fedemp/zonation-static-build/lib")
    target_link_libraries(zig4 zig4lib -Wl,-rpath -Wl,$ORIGIN/zigdso  -Wl,-as-needed -Wl,-dynamic-linker /lib64/ld-linux-x86-64.so.2 ${LOCAL_STATIC_LIBDIR}/libQtNetwork.a ${LOCAL_STATIC_LIBDIR}/libQtCore.a ${LOCAL_STATIC_LIBDIR}/libgdal.a    /usr/lib/x86_64-linux-gnu/libfftw3f.a /usr/lib/x86_64-linux-gnu/libfftw3f_omp.a /usr/lib/gcc/x86_64-linux-gnu/4.7/libgomp.a      /usr/lib/libboost_program_options.a /usr/lib/libboost_chrono.a /usr/lib/libboost_filesystem.a /usr/lib/libboost_system.a /usr/lib/libboost_thread.a           /usr/lib/libgeos_c.a /usr/lib/libgeos.a /usr/lib/x86_64-linux-gnu/libjasper.a     /usr/lib/libhdf5.a /usr/lib/libmfhdfalt.a /usr/lib/libdfalt.a /usr/lib/libgif.a /usr/lib/x86_64-linux-gnu/libjpeg.a /usr/lib/x86_64-linux-gnu/libpng.a /usr/lib/x86_64-linux-gnu/libz.a /usr/lib/x86_64-linux-gnu/libm.a  /usr/lib/x86_64-linux-gnu/librt.a -lpthread -ldl -Wl,-Bstatic -lm -Wl,-Bdynamic)  # -lpthread -ldl  #/usr/lib/gcc/x86_64-linux-gnu/4.7/libstdc++.a /usr/lib/x86_64-linux-gnu/libc.a       

  else()
    set(FFTW3_LIBRARIES fftw3f fftw3f_threads)
    target_link_libraries(zig4 zig4lib ${Boost_LIBRARIES} ${FFTW3_LIBRARIES} ${QT_LIBRARIES} ${GDAL_LIBRARIES} -lz)
  endif()
endif()

install(TARGETS zig4 DESTINATION bin)

# Here there should be some add_test...
