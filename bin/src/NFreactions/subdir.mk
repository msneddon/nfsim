################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFreactions/NFreactions.cpp 

OBJS += \
./src/NFreactions/NFreactions.o 

CPP_DEPS += \
./src/NFreactions/NFreactions.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFreactions/%.o: ../src/NFreactions/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


