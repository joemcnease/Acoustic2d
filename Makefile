#                                   _ ._  _ , _ ._
#                                 (_ ' ( `  )_  .__)
#                               ( (  (    )   `)  ) _)
#                              (__ (_   (_ . _) _) ,__)
#                                  `~~`\ ' . /`~~`
#                                  ,::: ;   ; :::,
#                                 ':::::::::::::::'
#  ___________________________jgs______/_ __ \___________________________________
# |                                                                              |
# | Makefile for building opengl/imgui application with glfw3/glad/glm/stbh/etc. |
# | Author: Joe McNease                                                          |
# | Date  : Tue 14 Dec 2021 06:29:19 PM CST                                      |
# |______________________________________________________________________________|


# Compiler details
CXX := g++ -std=c++17
CXXFLAGS := -g -Wall
LDFLAGS := -ldl -lGL -lpthread -lX11 -lXrandr -lglfw -lXi -lm

# Target file
TARGET := main

# Directories
SRCDIR := src
OBJDIR := obj

# All .cpp files should be in src to be built
SRCS := $(wildcard $(SRCDIR)/*.cpp)
OBJS := $(SRCS:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
OBJS += $(wildcard $(SRCDIR)/*.c)
DEPS := $(wildcard $(SRCDIR)/*.vs) $(wildcard $(SRCDIR)/*.fs)

# Recipes for disaster
all: $(TARGET)

# Build main by linking object files and ldflags to main
# Target depends on object files
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS)

# All src/.cpp files compiled to obj/.o files
# Object files depend on source files
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp $(DEPS) | $(OBJDIR)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Before the object files can be compiled we need to mkdir
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Clean up object files
.PHONY: clean
clean:
	rm -rf $(OBJDIR)

# Clean up all files made from build including $(TARGET)
.PHONY: clean-all
clean-all:
	rm -rf $(OBJDIR) $(TARGET)
