################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFutil/MTrand/mtrand.cpp \
../src/NFutil/MTrand/mttest.cpp 

OBJS += \
./src/NFutil/MTrand/mtrand.o \
./src/NFutil/MTrand/mttest.o 

CPP_DEPS += \
./src/NFutil/MTrand/mtrand.d \
./src/NFutil/MTrand/mttest.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFutil/MTrand/%.o: ../src/NFutil/MTrand/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include/c++/4.8.5 -O3 -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


