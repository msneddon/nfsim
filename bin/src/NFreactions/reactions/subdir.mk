################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFreactions/reactions/DORreaction.cpp \
../src/NFreactions/reactions/reaction.cpp 

OBJS += \
./src/NFreactions/reactions/DORreaction.o \
./src/NFreactions/reactions/reaction.o 

CPP_DEPS += \
./src/NFreactions/reactions/DORreaction.d \
./src/NFreactions/reactions/reaction.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFreactions/reactions/%.o: ../src/NFreactions/reactions/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


