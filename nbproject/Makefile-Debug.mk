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
	build/Debug/GNU-Linux-x86/load.o \
	build/Debug/GNU-Linux-x86/main.o \
	build/Debug/GNU-Linux-x86/efcm.o

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
.build-conf: ${BUILD_SUBPROJECTS} bin/sefcm

bin/sefcm: ${OBJECTFILES}
	${MKDIR} -p bin
	${LINK.cc} -o bin/sefcm ${OBJECTFILES} ${LDLIBSOPTIONS} 

build/Debug/GNU-Linux-x86/load.o: load.cpp 
	${MKDIR} -p build/Debug/GNU-Linux-x86
	$(COMPILE.cc) -g -o build/Debug/GNU-Linux-x86/load.o load.cpp

build/Debug/GNU-Linux-x86/main.o: main.cpp 
	${MKDIR} -p build/Debug/GNU-Linux-x86
	$(COMPILE.cc) -g -o build/Debug/GNU-Linux-x86/main.o main.cpp

build/Debug/GNU-Linux-x86/efcm.o: efcm.cpp 
	${MKDIR} -p build/Debug/GNU-Linux-x86
	$(COMPILE.cc) -g -o build/Debug/GNU-Linux-x86/efcm.o efcm.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf:
	${RM} -r build/Debug
	${RM} bin/sefcm

# Subprojects
.clean-subprojects:
