################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFoutput/NFoutput.cpp 

OBJS += \
./src/NFoutput/NFoutput.o 

CPP_DEPS += \
./src/NFoutput/NFoutput.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFoutput/%.o: ../src/NFoutput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


