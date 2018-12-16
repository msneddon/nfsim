################################################################################
# Automatically-generated file. Do not edit!
################################################################################

-include ../makefile.init

RM := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include src/nauty24/subdir.mk
-include src/NFutil/MTrand/subdir.mk
-include src/NFutil/subdir.mk
-include src/NFtest/transcription/subdir.mk
-include src/NFtest/tlbr/subdir.mk
-include src/NFtest/simple_system/subdir.mk
-include src/NFtest/agentcell/cell/subdir.mk
-include src/NFtest/agentcell/subdir.mk
-include src/NFscheduler/subdir.mk
-include src/NFreactions/transformations/subdir.mk
-include src/NFreactions/reactions/subdir.mk
-include src/NFreactions/reactantLists/subdir.mk
-include src/NFreactions/mappings/subdir.mk
-include src/NFreactions/subdir.mk
-include src/NFoutput/subdir.mk
-include src/NFinput/TinyXML/subdir.mk
-include src/NFinput/subdir.mk
-include src/NFfunction/muParser/subdir.mk
-include src/NFfunction/subdir.mk
-include src/NFcore/reactionSelector/subdir.mk
-include src/NFcore/moleculeLists/subdir.mk
-include src/NFcore/subdir.mk
-include src/subdir.mk
-include subdir.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: NFsim

# Tool invocations
NFsim: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o "NFsim" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) NFsim
	-@echo ' '

.PHONY: all clean dependents

-include ../makefile.targets
