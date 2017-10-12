################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFtest/transcription/transcription.cpp 

OBJS += \
./src/NFtest/transcription/transcription.o 

CPP_DEPS += \
./src/NFtest/transcription/transcription.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFtest/transcription/%.o: ../src/NFtest/transcription/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


