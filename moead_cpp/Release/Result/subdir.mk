################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Result/IGD.cpp \
../Result/cec09filter.cpp 

OBJS += \
./Result/IGD.o \
./Result/cec09filter.o 

CPP_DEPS += \
./Result/IGD.d \
./Result/cec09filter.d 


# Each subdirectory must supply rules for building sources it contributes
Result/%.o: ../Result/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


