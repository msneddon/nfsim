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
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


