PROG_NAME   = $(notdir $(CURDIR))
DIST_DIR    = ./dist

SRC_DIR     = ./src
BUILD_DIR   = ./build
INCLUDE_DIR = ./include

ifeq ($(OS), Windows_NT)
	LIB_DIR = ./lib/
else
	LIB_DIR = ./liblinux/
endif

CFLAGS      = -Wall -O3 -march=native -fopenmp
LFLAGS      = -static -fopenmp
LIBS        = -L$(LIB_DIR) -lraylib -lm

ifeq ($(OS), Windows_NT)
	LIBS += -lopengl32 -lgdi32 -lwinmm
endif

CC          = gcc
LD          = gcc

rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))

SRC_LIST = $(call rwildcard, $(SRC_DIR), *.c)
OBJ_LIST = $(addsuffix .o, $(basename $(patsubst %, $(BUILD_DIR)/%, $(SRC_LIST:$(SRC_DIR)/%=%))))

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) -c $(CFLAGS) -I $(INCLUDE_DIR) -o $@ $<
	
#compile and link
default: $(PROG_NAME)

#compile without linking
compile: $(OBJ_LIST)

#compile and link
$(PROG_NAME): compile
	$(LD) $(LFLAGS) $(OBJ_LIST) $(LIBS) -o $(DIST_DIR)/$@

clean:
	rm -f $(DIST_DIR)/$(PROG_NAME) $(BUILD_DIR)/*.o