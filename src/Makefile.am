ACLOCAL_AMFLAGS = -I m4
AUTOMAKE_OPTIONS = subdir-objects

bin_PROGRAMS=coverageDropping coverageDropping_debug


coverageDropping_SOURCES=coverageDropping.cpp \
	data_structures/Mask.cpp \
	data_structures/Contig.cpp \
	errors/Data_Exception.cpp \
	errors/File_Not_Found.cpp \
	errors/Incorrect_Format.cpp \
	io/Auto_Unzip.cpp \
	io/Fasta.cpp \
	io/read_functions.cpp \
	data_structures/Reads.cpp 

coverageDropping_debug_SOURCES=$(coverageDropping_SOURCES)

SAMTOOLS=samtools/
INCLUDES = -I$(SAMTOOLS) 

coverageDropping_CXXFLAGS= -O3 $(BOOST_CPPFLAGS)
coverageDropping_debug_CXXFLAGS= -O0 -g3 $(BOOST_CPPFLAGS)
coverageDropping_CFLAGS= -O3 
coverageDropping_debug_CFLAGS= -O0 -g3 


coverageDropping_LDFLAGS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS)  $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
coverageDropping_LIBS = -L$(SAMTOOLS) -lbam $(BOOST_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_IOSTREAMS_LIBS)  -lz -lbz2
coverageDropping_debug_LDFLAGS = $(coverageDropping_LDFLAGS)
coverageDropping_debug_LIBS = $(coverageDropping_LIBS)

