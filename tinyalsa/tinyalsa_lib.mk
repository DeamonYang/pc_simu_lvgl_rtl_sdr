#CSRCS += $(shell find -L src -name "*.c")

CSRCS += src/mixer_hw.c
CSRCS += src/mixer_plugin.c
CSRCS += src/pcm_hw.c
CSRCS += src/pcm_plugin.c
CSRCS += src/limits.c
CSRCS += src/pcm.c
CSRCS += src/snd_card_plugin.c
CSRCS += src/mixer.c



CFLAGS += "-I$(LVGL_DIR)/tinyalsa/src"
