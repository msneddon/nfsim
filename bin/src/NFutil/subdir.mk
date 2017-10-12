################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFutil/conversion.cpp \
../src/NFutil/random.cpp \
../src/NFutil/stringOperations.cpp 

OBJS += \
./src/NFutil/conversion.o \
./src/NFutil/random.o \
./src/NFutil/stringOperations.o 

CPP_DEPS += \
./src/NFutil/conversion.d \
./src/NFutil/random.d \
./src/NFutil/stringOperations.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFutil/%.o: ../src/NFutil/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


