OBROOT = D:\openbabel-2.3.2
DESROOT = D:\cpplibrary\descriptors
CXX = cl.exe
CFLAGS = /I $(OBROOT)\include\openbabel-2.0 /I $(DESROOT)/include /W3 /O2 /MD /D NDEBUG /EHsc /D WIN32 /D USING_DYNAMIC_LIBS
LDFLAGS = /INCREMENTAL:NO /SUBSYSTEM:CONSOLE /LIBPATH:$(OBROOT)/bin /LIBPATH:$(DESROOT)/bin openbabel-2.lib descriptors.lib

all: calcDescriptors.exe


main.obj: main.cpp
	$(CXX) $(CFLAGS) $< /c

calcDescriptors.exe: main.obj
	link.exe main.obj $(LDFLAGS) /OUT:calcDescriptors.exe

