################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFtest/agentcell/agentcell.cpp 

OBJS += \
./src/NFtest/agentcell/agentcell.o 

CPP_DEPS += \
./src/NFtest/agentcell/agentcell.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFtest/agentcell/%.o: ../src/NFtest/agentcell/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


