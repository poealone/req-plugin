# ReQ — PocketDAW Parametric EQ
# Usage: make [linux|anbernic|windows|all]

SRCS      = req.c
CFLAGS    = $(shell sdl2-config --cflags 2>/dev/null)
LDFLAGS   = -lm $(shell sdl2-config --libs 2>/dev/null)

linux: $(SRCS)
	gcc -shared -fPIC -O2 $(CFLAGS) -o req.so $(SRCS) $(LDFLAGS)
	@echo "Built: req.so (Linux x64)"

anbernic: $(SRCS)
	aarch64-none-linux-gnu-gcc -shared -fPIC -O2 -I/usr/include/SDL2 -o req.so $(SRCS) -lm -lSDL2
	@echo "Built: req.so (Anbernic aarch64)"

windows: $(SRCS)
	x86_64-w64-mingw32-gcc -shared -I/usr/x86_64-w64-mingw32/include/SDL2 -o req.dll $(SRCS) -lm -lSDL2
	@echo "Built: req.dll (Windows x64)"

all: linux

.PHONY: linux anbernic windows all
