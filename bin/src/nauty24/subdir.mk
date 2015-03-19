################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/nauty24/nausparse.c \
../src/nauty24/nautil.c \
../src/nauty24/nauty.c 

OBJS += \
./src/nauty24/nausparse.o \
./src/nauty24/nautil.o \
./src/nauty24/nauty.o 

C_DEPS += \
./src/nauty24/nausparse.d \
./src/nauty24/nautil.d \
./src/nauty24/nauty.d 


# Each subdirectory must supply rules for building sources it contributes
src/nauty24/%.o: ../src/nauty24/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


