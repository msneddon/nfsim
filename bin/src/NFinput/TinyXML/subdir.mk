################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/NFinput/TinyXML/tinystr.cpp \
../src/NFinput/TinyXML/tinyxml.cpp \
../src/NFinput/TinyXML/tinyxmlerror.cpp \
../src/NFinput/TinyXML/tinyxmlparser.cpp 

OBJS += \
./src/NFinput/TinyXML/tinystr.o \
./src/NFinput/TinyXML/tinyxml.o \
./src/NFinput/TinyXML/tinyxmlerror.o \
./src/NFinput/TinyXML/tinyxmlparser.o 

CPP_DEPS += \
./src/NFinput/TinyXML/tinystr.d \
./src/NFinput/TinyXML/tinyxml.d \
./src/NFinput/TinyXML/tinyxmlerror.d \
./src/NFinput/TinyXML/tinyxmlparser.d 


# Each subdirectory must supply rules for building sources it contributes
src/NFinput/TinyXML/%.o: ../src/NFinput/TinyXML/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


