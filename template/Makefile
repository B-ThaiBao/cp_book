PROBLEM_NAME := problem_name
LANG := lang
DEBUG := true

SRC := $(PROBLEM_NAME).$(LANG)
DEBUG_FLAG := $(if $(filter true 1 yes on,$(DEBUG)),-d,)

.PHONY: all run test build clean

# This is for nvim cmd
all:
	@if ls *.in >/dev/null 2>&1; then maketest $(DEBUG_FLAG) -ioc $(SRC); else makerun $(DEBUG_FLAG) -c $(SRC); fi

# This is for normal terminal
run:
	@makerun $(DEBUG_FLAG) $(SRC)

test:
	@maketest $(DEBUG_FLAG) -io $(SRC)

build:
	@makebuild $(DEBUG_FLAG) $(SRC)

clean:
	@makeclean
