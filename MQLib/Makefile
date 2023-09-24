CFLAGS = -Iinclude -std=c++0x -O2 -Wall
LFLAGS = -lm
SHAREFLAGS = -shared -fPIC
BUILDDIR = .build
SRCDIR = src
EXECUTABLE = bin/MQLib
STATIC = bin/MQLib.a
SRCS = $(shell find $(SRCDIR) -name "*.cpp")
OBJS = $(shell echo "$(SRCS)" | sed -e "s/ $(SRCDIR)/ $(BUILDDIR)/g" -e "s/^$(SRCDIR)/$(BUILDDIR)/g" -e "s/\.cpp/.o/g")
DEPS = $(shell echo "$(OBJS)" | sed -e "s/\.o/.P/g")

all: $(EXECUTABLE) $(STATIC)

$(EXECUTABLE): $(OBJS)
	g++ $(LFLAGS) -o $(EXECUTABLE) $(OBJS)

$(STATIC): $(OBJS)
	ar -cq $(STATIC) $(OBJS)

### Conversion from .d to .P from http://mad-scientist.net/make/autodep.html
$(BUILDDIR)/%.o : $(SRCDIR)/%.cpp
	g++ -c -MD $(CFLAGS) -o $@ $<
	@cp $(BUILDDIR)/$(*).d $(BUILDDIR)/$(*).P; \
	sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' -e '/^$$/ d' -e 's/$$/ :/' < $(BUILDDIR)/$(*).d >> $(BUILDDIR)/$(*).P; \
	rm -f $(BUILDDIR)/$(*).d

clean:
	@rm -f $(OBJS) $(DEPS) $(EXECUTABLE) $(STATIC)
	@rm -f `find . -name "*~"`
	@rm -f `find . -name ".DS_Store"`
	@rm -f `find . -name "*\.aux"`
	@rm -f `find . -name "*\.log"`
	@rm -f `find . -name "*\.synctex\.gz"`
	@rm -f `find . -name "*\.Rapp\.history"`

-include $(DEPS)
