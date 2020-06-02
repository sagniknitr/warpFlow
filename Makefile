# Use Bash sdf
SHELL = /bin/sh
 
# Functions
find_includes_in_dir = $(shell find $(1) -name "*.hpp" | sed 's|/[^/]*$$||' | sort -u)
 
# ---------------------------------------------------------------------
# Toolchain Configuration
# ---------------------------------------------------------------------
C_STANDARD              := -std=gnu11
CXX_STANDARD            := -std=gnu++14
 
# -----------------------------------------------------------------------------------------------------------------
# Defined Symbols
# -----------------------------------------------------------------------------------------------------------------
DEFS                    := 
 
# ---------------------------------------------------------------------------------------------------------------------------------------
# Compiler & Linker Flags
# ---------------------------------------------------------------------------------------------------------------------------------------
# Flags sent to all tools in the Toolchain 
TOOLCHAIN_SETTINGS      := -fmessage-length=0
 
# C Compiler -- Warnings 
CFLAGS                  += $(TOOLCHAIN_SETTINGS) $(DEFS) $(addprefix -I, $(INC_DIRS))
CFLAGS                  += -Wall
CFLAGS                  += -Wextra
CFLAGS                  += -Wfatal-errors
CFLAGS                  += -Wpacked
CFLAGS                  += -Winline
CFLAGS                  += -Wfloat-equal
CFLAGS                  += -Wconversion
CFLAGS                  += -Wpointer-arith
CFLAGS                  += -Wdisabled-optimization
CFLAGS                  += -Wno-unused-parameter
 
# C++ Compiler -- Required & Optimization Flags
CXXFLAGS                += $(CFLAGS)
 
# C++ -- Warnings
CXXFLAGS                += -Weffc++
CXXFLAGS                += -Wfloat-equal
CXXFLAGS                += -Wsign-promo
CXXFLAGS                += -Wmissing-declarations 
CXXFLAGS                += -Woverloaded-virtual
CXXFLAGS                += -Wmissing-format-attribute
CXXFLAGS                += -Wold-style-cast
CXXFLAGS                += -Wshadow
CXXFLAGS                += -Wctor-dtor-privacy
 
# Linker
LDFLAGS                 += $(TOOLCHAIN_SETTINGS) $(DEFS)
 
# -------------------------------------------------------------
# Build Type Modifiers
# -------------------------------------------------------------
# Debug
DEFS_DEBUG              := -DDEBUG
CFLAGS_DEBUG            := -ggdb -g3 -Og
 
# Release
CFLAGS_RELEASE          := -O3
 
#########################################################################################################################################
# RULE DEFINITIONS -- This section is generic
#########################################################################################################################################
 
# =======================================================================================================================================
# Build Configuration Rule 
# - Generate build config using Product Root Directory ($1), Build Type ("Debug" or "Release") ($2)
# =======================================================================================================================================
define CONFIG_RULE
BUILD_DIR               := $1/Build/$2
OBJ_DIR                 := $$(BUILD_DIR)/obj
INC_DIRS                := $$(call find_includes_in_dir, $$(SRC_DIRS))
HEADERS                 := $$(foreach dir, $$(SRC_DIRS), $$(shell find $$(dir) -name "*.hpp"))
ASM_SRC                 := $$(foreach dir, $$(SRC_DIRS), $$(shell find $$(dir) -name "*.s"))
C_SRC                   := $$(foreach dir, $$(SRC_DIRS), $$(shell find $$(dir) -name "*.c"))
CXX_SRC                 := $$(foreach dir, $$(SRC_DIRS), $$(shell find $$(dir) -name "*.cpp"))
OBJECTS                 := $$(addprefix $$(OBJ_DIR)/, $$(C_SRC:.c=.o) $$(CXX_SRC:.cpp=.o) $$(ASM_SRC:.s=.o))
LDSCRIPTS               := $$(addprefix -T, $$(foreach dir, $$(SRC_DIRS), $$(shell find $$(dir) -name "*.ld")))
DIRS                    := $$(BUILD_DIR) $$(sort $$(dir $$(OBJECTS)))
AUTODEPS                := $$(OBJECTS:.o=.d)
 
 
ifeq ($2, Release)
    DEFS    += $$(DEFS_RELEASE)
    CFLAGS  += $$(CFLAGS_RELEASE)
    LDFLAGS += $$(LDFLAGS_RELEASE)
else 
    DEFS    += $$(DEFS_DEBUG)
    CFLAGS  += $$(CFLAGS_DEBUG)
    LDFLAGS += $$(LDFLAGS_DEBUG)
endif
 
endef 
# =======================================================================================================================================
# End CONFIG_RULE
# =======================================================================================================================================
 
 
# =======================================================================================================================================
# Build Target Rule 
# - Generate build config using Product Name ($1), Product Root Directory ($2), Build Type ("Debug" or "Release") ($3)
# =======================================================================================================================================
define BUILD_TARGET_RULE
$(eval $(call CONFIG_RULE,$2,$3))
 
all : $$(BUILD_DIR)/$1
 
# Tool Invocations
$$(BUILD_DIR)/$1 : $$(OBJECTS) | $$(BUILD_DIR)
    @echo ' '
    @echo 'Building $$(@)'
    @echo 'Invoking: C++ Linker'
    $$(CXX) $$(LDFLAGS) $$(LDSCRIPTS) -o $$(@) $$(OBJECTS)
    @echo 'Finished building: $$@'
    @echo ' '
 
$$(OBJECTS) : | $$(DIRS)
 
$$(DIRS) : 
    @echo Creating $$(@)
    @mkdir -p $$(@)
 
$$(OBJ_DIR)/%.o : %.c
    @echo Compiling $$(<F)
    @$$(CC) $$(C_STANDARD) $$(CFLAGS) -c -MMD -MP $$< -o $$(@)
 
$$(OBJ_DIR)/%.o : %.cpp
    @echo Compiling $$(<F)
    @$$(CXX) $$(CXX_STANDARD) $$(CXXFLAGS) -c -MMD -MP $$< -o $$(@)
 
$$(OBJ_DIR)/%.o : %.s
    @echo Assembling $$(<F)
    @$$(AS) $$(ASFLAGS) $$< -o $$(@)
 
clean :
    @rm -rf $$(PRODUCT_DIR)/Build
 
.PHONY : clean all
 
# include by auto dependencies
-include $$(AUTODEPS)
 
endef
# =======================================================================================================================================
# End BUILD_TARGET_RULE
# =======================================================================================================================================
#########################################################################################################################################
#########################################################################################################################################
 
# Build Type
ifeq ($(build), Debug)
    BUILD_TYPE := Debug
else
    BUILD_TYPE := Release
endif
 
 
# Defaults
PRODUCT ?= OpticalFlow_Warping
PRODUCT_DIR ?= .
BUILD_TYPE ?= Debug
SRC_DIRS ?= ./Source/cpp
 
# Evaluate Rules Defined Above
$(eval $(call BUILD_TARGET_RULE,$(PRODUCT),$(PRODUCT_DIR),$(BUILD_TYPE)))
