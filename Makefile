CC      = clang++
CFLAGS  = -g -std=c++11 -O2
LDFLAGS = 
LIBS    =  /usr/local/Cellar/boost/1.67.0_1/lib
INC     = ./hpp /usr/local/Cellar/boost/1.67.0_1/include
INC_PARAMS=$(foreach d, $(INC), -I$d)
SRC_DIR = ./src
OBJ_DIR = ./build
SOURCES = $(shell ls $(SRC_DIR)/*.cpp) 
OBJS    = $(subst $(SRC_DIR),$(OBJ_DIR), $(SOURCES:.cpp=.o))
TARGET  = a.out
DEPENDS = $(OBJS:.o=.d)
all: $(TARGET)
$(TARGET): $(OBJS) $(LIBS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp 
	@if [ ! -d $(OBJ_DIR) ]; \
		then echo "mkdir -p $(OBJ_DIR)"; mkdir -p $(OBJ_DIR); \
		fi
	$(CC) $(CFLAGS) $(INC_PARAMS) -o $@ -c $< 
clean:
	$(RM) $(OBJS) $(TARGET) $(DEPENDS)
-include $(DEPENDS)
.PHONY: all clean
