################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFinput/NFinput.cpp \
../src/NFinput/commandLineParser.cpp \
../src/NFinput/parseFuncXML.cpp \
../src/NFinput/parseSymRxns.cpp \
../src/NFinput/rnfRunner.cpp \
../src/NFinput/walk.cpp 

OBJS += \
./src/NFinput/NFinput.o \
./src/NFinput/commandLineParser.o \
./src/NFinput/parseFuncXML.o \
./src/NFinput/parseSymRxns.o \
./src/NFinput/rnfRunner.o \
./src/NFinput/walk.o 

CPP_DEPS += \
./src/NFinput/NFinput.d \
./src/NFinput/commandLineParser.d \
./src/NFinput/parseFuncXML.d \
./src/NFinput/parseSymRxns.d \
./src/NFinput/rnfRunner.d \
./src/NFinput/walk.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFinput/%.o: ../src/NFinput/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


