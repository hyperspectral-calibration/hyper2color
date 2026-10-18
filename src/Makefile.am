bin_PROGRAMS =	hyper2color
AM_CFLAGS = @OPENMP_CFLAGS@

LIBS = @LIBS@ -ltiff -lz -ljpeg
hyper2color_SOURCES = \
			hyspex.h \
			hyspex.c \
			CIE-tristimulus.c \
			CIE-A.c \
			CIE-D65.c \
			icc.c \
			color.h \
			color.c \
			hyper2color.c
