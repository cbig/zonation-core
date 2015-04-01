cmake_minimum_required(VERSION 2.8)
# Intended usage:
# For normal (compact) vmat structure (COMPACT_VMAT_LOOPUP defined by default):
#      cmake ../path/to/zonation/src/
# win: cmake -G"MSYS Makefiles" ../path/to/zonation/src/
# For old style / traditional / non-compact vmat structure (COMPACT_VMAT_FLOATP definition required):
#      cmake -DCOMPACT_VMAT_FLOATP=1 ../path/to/zonation/src/
# win: cmake -DCOMPACT_VMAT_FLOATP=1 -G"MSYS Makefiles" ../path/to/zonation/src/
#      for 32 bits compilation add -DMINGW_IS_W32=1
project(zonation)

set(zonation_VERSION_MAJOR 4)
set(zonation_VERSION_MINOR 0)
set(zonation_VERSION_PATCH 0)

if (NOT COMPACT_VMAT_LOOKUP)
   set(COMPACT_VMAT_LOOKUP 1)
endif()

if (COMPACT_VMAT_LOOKUP)
  set(zonation_VERSION_PATCH "${zonation_VERSION_PATCH}_compact")
else()
  set(zonation_VERSION_PATCH "${zonation_VERSION_PATCH}_traditional")
endif()

# For specificities of w64 as opposed to standard 32 bits MinGW
# For 64 bits build: comment it out or define it with a -D command line parameter
if (WIN32 AND NOT MINGW_IS_W32)
   set(MINGW_IS_W64 1)
endif()
# If there're problems with make missing in some MinGW distribs
# set(CMAKE_MAKE_PROGRAM /c/local/MinGW/bin/mingw32-make.exe)

#include_directories(/local/include/)

if (COMPILE_STATIC)
  message("*************** Note: compiling static! ***************")
  #set(CMAKE_EXE_LINKER_FLAGS "-static")    # Problems with Qt: getpwnam_r, getaddrinfo, gethostbyaddr, dlopen, etc
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} "-static-libgcc -static-libstdc++")
endif()


# if compiling static you'll need this
# boost thread doc does not even mention this!
if (MINGW)
   message("***** Note: MINGW.*****")
   set(LOCAL_LIBDIR /mingw/local/lib)
   set(CMAKE_EXE_LINKER_FLAGS "-static") 

   add_definitions(-DBOOST_THREAD_USE_LIB)
   # for dynamic compililng it would be -DBOOST_THREAD_USE_DLL

   set(CXX_OPT_OPT "-O2 -finline-functions -march=core2 -mfpmath=sse -mieee-fp -funroll-loops -pipe -mmmx -msse -msse2 -msse3 -m3dnow")  # More risk: -mlwp -mabm -mbmi -mtbm    / i3/i5/i7: avx!
   # Note: avoid -ffast-math, it breaks NaN things.
   if (MINGW_IS_W64)
      message("***** Note: configuring for MinGW-w64 (64 bits + 32 bits (cross)). Make sure that's right.*****")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_OPT_OPT}")  
      set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} "-static-libgcc")
   else()
      message("***** Note: configuring for MinGW (32 bits). Make sure that's right.*****")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_OPT_OPT}")
      set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} "-static-libgcc -static-libstdc++")
   endif()

   set(Boost_LIBRARIES "/local/lib/libboost_system-mgw46-mt-1_48.a  /local/lib/libboost_filesystem-mgw46-mt-1_48.a /local/lib/libboost_program_options-mgw46-mt-1_48.a /local/lib/libboost_thread-mgw46-mt-1_48.a /local/lib/libboost_signals-mgw46-mt-1_48.a") 
endif()

configure_file("${PROJECT_SOURCE_DIR}/config.h.in" "${PROJECT_BINARY_DIR}/cmake_generated/config.h")
include_directories("${PROJECT_BINARY_DIR}/cmake_generated/")

# Checks for dependencies inside the subdirs...
add_subdirectory(zig4lib)
add_subdirectory(zig4)
