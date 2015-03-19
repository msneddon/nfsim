################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFfunction/compositeFunction.cpp \
../src/NFfunction/funcParser.cpp \
../src/NFfunction/function.cpp \
../src/NFfunction/localFunction.cpp 

OBJS += \
./src/NFfunction/compositeFunction.o \
./src/NFfunction/funcParser.o \
./src/NFfunction/function.o \
./src/NFfunction/localFunction.o 

CPP_DEPS += \
./src/NFfunction/compositeFunction.d \
./src/NFfunction/funcParser.d \
./src/NFfunction/function.d \
./src/NFfunction/localFunction.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFfunction/%.o: ../src/NFfunction/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


