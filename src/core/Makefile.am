lib_LIBRARIES = libcore.a

libcore_a_SOURCES = \
  electric_eqn.f90 \
  magnetic_eqn.f90 \
  calc_dt.f90 \
  e2exyz.f90 \
  array_add.f90 \
  calc_poynting.f90 \
  boundary.cxx
libcore_a_FCFLAGS = 
libcore_a_CFLAGS   = -I ../include
libcore_a_CXXFLAGS = -I ../include

check_PROGRAMS    = core_e2m
core_e2m_SOURCES  = test/e2m.cxx
core_e2m_LDADD    = libcore.a
core_e2m_CFLAGS   = -I ../include
core_e2m_CXXFLAGS = -I ../include

TESTS = core_e2m
dist_noinst_DATA = $(TEST)

