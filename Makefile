#
# Makefile
#
$(warning Using Make to build this project is deprecated, please switch to CMake)
CC ?= gcc
LVGL_DIR_NAME ?= lvgl
LVGL_DIR ?= ${shell pwd}

WARNINGS ?= -Wall -Wextra \
						-Wshadow -Wundef -Wmaybe-uninitialized -Wmissing-prototypes -Wno-discarded-qualifiers \
						-Wno-unused-function -Wno-error=strict-prototypes -Wpointer-arith -fno-strict-aliasing -Wno-error=cpp -Wuninitialized \
						-Wno-unused-parameter -Wno-missing-field-initializers -Wno-format-nonliteral -Wno-cast-qual -Wunreachable-code -Wno-switch-default  \
					  -Wreturn-type -Wmultichar -Wformat-security -Wno-ignored-qualifiers -Wno-error=pedantic -Wno-sign-compare -Wno-error=missing-prototypes -Wdouble-promotion -Wclobbered -Wdeprecated  \
						-Wempty-body  -Wstack-usage=2048 \
            -Wtype-limits -Wsizeof-pointer-memaccess -Wpointer-arith
            
CFLAGS ?= -O3 -I$(LVGL_DIR)/ $(WARNINGS) -I ./rtl_sdr_app_work/include/ `pkg-config --cflags wayland-client` `pkg-config --cflags xkbcommon`
LDFLAGS ?= -lSDL2 -lpthread -lusb-1.0 -lm -lfftw3 `pkg-config --libs wayland-client` `pkg-config --libs xkbcommon` -lpthread
BIN = demo


#Collect the files to compile
MAINSRC = ./main.c

include $(LVGL_DIR)/lvgl/lvgl.mk
include $(LVGL_DIR)/lv_drivers/lv_drivers.mk
include $(LVGL_DIR)/lv_demos/lv_demo.mk

CSRCS +=$(LVGL_DIR)/mouse_cursor_icon.c
CSRCS += ./rtl_sdr_app_work/lib/librtlsdr.c
CSRCS += ./rtl_sdr_app_work/lib/tuner_e4k.c
CSRCS += ./rtl_sdr_app_work/lib/tuner_fc0012.c
CSRCS += ./rtl_sdr_app_work/lib/tuner_fc0013.c
CSRCS += ./rtl_sdr_app_work/lib/tuner_fc2580.c
CSRCS += ./rtl_sdr_app_work/lib/tuner_r82xx.c
CSRCS += ./rtl_sdr_app_work/convenience/convenience.c

RTLSDR_SOURCE=$(wildcard ./rtl_sdr_app_work/lib/*.c  ./rtl_sdr_app_work/convenience/*.c)

OBJEXT ?= .o

AOBJS = $(ASRCS:.S=$(OBJEXT))
COBJS = $(CSRCS:.c=$(OBJEXT))

MAINOBJ = $(MAINSRC:.c=$(OBJEXT))

SRCS = $(ASRCS) $(CSRCS) $(MAINSRC)
OBJS = $(AOBJS) $(COBJS)


## MAINOBJ -> OBJFILES

all: default

%.o: %.c
	@$(CC)  $(CFLAGS) -c $< -o $@
	@echo "CC $<"
    
default: $(AOBJS) $(COBJS) $(MAINOBJ) 
	$(CC) -o $(BIN) $(MAINOBJ) $(AOBJS) $(COBJS) $(LDFLAGS)

clean: 
	rm -f $(BIN) $(AOBJS) $(COBJS) $(MAINOBJ)

