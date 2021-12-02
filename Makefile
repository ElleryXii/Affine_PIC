# This is for GNU make; other versions of make may not run correctly.

MAIN_PROGRAM = flip2d
SRC = grid.cpp particles.cpp main.cpp
MAIN_WITH_VIEWER = flip2dv
SRC_WITH_VIEWER = grid.cpp particles.cpp mainwithviewer.cpp viewflip2d/gluvi.cpp

include Makefile.defs

# object files
RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC)))
DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC)))


# object files
RELEASE_OBJ_V = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC_WITH_VIEWER)))
DEBUG_OBJ_V = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC_WITH_VIEWER)))

# how to make the main target (debug mode, the default)
$(MAIN_PROGRAM): $(DEBUG_OBJ)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS)

# how to make the main target (release mode)
$(MAIN_PROGRAM)_release: $(RELEASE_OBJ)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS)


# how to make the main target (debug mode, the default)
$(MAIN_WITH_VIEWER): $(DEBUG_OBJ_V)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS_V) $(GL_FLAGS)

# how to make the main target (release mode)
$(MAIN_WITH_VIEWER)_release: $(RELEASE_OBJ_V)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS_V) $(GL_FLAGS)

.PHONY: release
release: $(MAIN_PROGRAM)_release

.PHONY: debug
debug: $(MAIN_PROGRAM)

.PHONY: release_v
release_v: $(MAIN_WITH_VIEWER)_release

.PHONY: debug_v
debug_v: $(MAIN_WITH_VIEWER)

# how to compile each file
.SUFFIXES:
obj/%.o:
	$(CC) -c $(RELEASE_FLAGS) -o $@ $<
obj_debug/%.o:
	$(CC) -c $(DEBUG_FLAGS) -o $@ $<

# cleaning up
.PHONY: clean
clean:
	-rm -f obj/*.o $(MAIN_PROGRAM) obj_debug/*.o $(MAIN_PROGRAM)_release *core viewer $(MAIN_WITH_VIEWER) $(MAIN_WITH_VIEWER)_release

# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir obj
	-rm -f obj/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) >> obj/depend;)
	$(foreach srcfile,$(SRC_WITH_VIEWER),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) >> obj/depend;)
	-mkdir obj_debug
	-rm -f obj_debug/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) >> obj_debug/depend;)
	$(foreach srcfile,$(SRC_WITH_VIEWER),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) >> obj_debug/depend;)

-include obj/depend
-include obj_debug/depend

