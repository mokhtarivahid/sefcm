#
# Gererated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add custumized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=g77

# Include project Makefile
include Makefile

# Object Files
OBJECTFILES= \
	build/Release/GNU-Linux-x86/load.o \
	build/Release/GNU-Linux-x86/main.o \
	build/Release/GNU-Linux-x86/efcm.o

# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS} dist/Release/GNU-Linux-x86/sefcm

dist/Release/GNU-Linux-x86/sefcm: ${OBJECTFILES}
	${MKDIR} -p dist/Release/GNU-Linux-x86
	${LINK.cc} -o dist/Release/GNU-Linux-x86/sefcm ${OBJECTFILES} ${LDLIBSOPTIONS} 

build/Release/GNU-Linux-x86/load.o: load.cpp 
	${MKDIR} -p build/Release/GNU-Linux-x86
	$(COMPILE.cc) -O2 -o build/Release/GNU-Linux-x86/load.o load.cpp

build/Release/GNU-Linux-x86/main.o: main.cpp 
	${MKDIR} -p build/Release/GNU-Linux-x86
	$(COMPILE.cc) -O2 -o build/Release/GNU-Linux-x86/main.o main.cpp

build/Release/GNU-Linux-x86/efcm.o: efcm.cpp 
	${MKDIR} -p build/Release/GNU-Linux-x86
	$(COMPILE.cc) -O2 -o build/Release/GNU-Linux-x86/efcm.o efcm.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Release
	${RM} dist/Release/GNU-Linux-x86/sefcm

# Subprojects
.clean-subprojects:
