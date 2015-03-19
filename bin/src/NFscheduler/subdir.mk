################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFscheduler/NFstream.cpp \
../src/NFscheduler/Scheduler.cpp 

OBJS += \
./src/NFscheduler/NFstream.o \
./src/NFscheduler/Scheduler.o 

CPP_DEPS += \
./src/NFscheduler/NFstream.d \
./src/NFscheduler/Scheduler.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFscheduler/%.o: ../src/NFscheduler/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


