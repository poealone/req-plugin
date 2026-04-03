# ReQ — PocketDAW Parametric EQ
# Usage: make [linux|anbernic|windows|all]

SRCS    = req.c
LDFLAGS = -lm

linux: $(SRCS)
	gcc -shared -fPIC -O2 -o req.so $(SRCS) $(LDFLAGS)
	@echo "Built: req.so (Linux x64)"

anbernic: $(SRCS)
	aarch64-none-linux-gnu-gcc -shared -fPIC -O2 -o req.so $(SRCS) $(LDFLAGS)
	@echo "Built: req.so (Anbernic aarch64)"

windows: $(SRCS)
	x86_64-w64-mingw32-gcc -shared -o req.dll $(SRCS) $(LDFLAGS)
	@echo "Built: req.dll (Windows x64)"

all: linux

.PHONY: linux anbernic windows all
