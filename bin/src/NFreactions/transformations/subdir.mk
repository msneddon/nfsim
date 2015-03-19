################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFreactions/transformations/moleculeCreator.cpp \
../src/NFreactions/transformations/speciesCreator.cpp \
../src/NFreactions/transformations/transformation.cpp \
../src/NFreactions/transformations/transformationSet.cpp 

OBJS += \
./src/NFreactions/transformations/moleculeCreator.o \
./src/NFreactions/transformations/speciesCreator.o \
./src/NFreactions/transformations/transformation.o \
./src/NFreactions/transformations/transformationSet.o 

CPP_DEPS += \
./src/NFreactions/transformations/moleculeCreator.d \
./src/NFreactions/transformations/speciesCreator.d \
./src/NFreactions/transformations/transformation.d \
./src/NFreactions/transformations/transformationSet.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFreactions/transformations/%.o: ../src/NFreactions/transformations/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


