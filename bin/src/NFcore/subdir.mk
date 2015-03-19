################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFcore/complex.cpp \
../src/NFcore/complexList.cpp \
../src/NFcore/molecule.cpp \
../src/NFcore/moleculeType.cpp \
../src/NFcore/observable.cpp \
../src/NFcore/reactionClass.cpp \
../src/NFcore/system.cpp \
../src/NFcore/templateMolecule.cpp 

OBJS += \
./src/NFcore/complex.o \
./src/NFcore/complexList.o \
./src/NFcore/molecule.o \
./src/NFcore/moleculeType.o \
./src/NFcore/observable.o \
./src/NFcore/reactionClass.o \
./src/NFcore/system.o \
./src/NFcore/templateMolecule.o 

CPP_DEPS += \
./src/NFcore/complex.d \
./src/NFcore/complexList.d \
./src/NFcore/molecule.d \
./src/NFcore/moleculeType.d \
./src/NFcore/observable.d \
./src/NFcore/reactionClass.d \
./src/NFcore/system.d \
./src/NFcore/templateMolecule.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFcore/%.o: ../src/NFcore/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


