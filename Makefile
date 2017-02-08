CC := gcc 

CFLAGS := -g -Wall -O2

CXXFLAGS := -std=c++0x -O3  -fno-exceptions   -Wall -Wextra 

OBJS = compress.o decompress.o lwfqzip.o 
OBJS1 = lwmapping.o assistant.o genProcess.o map_mulPreKey.o
OBJS2 = FQZip.cpp


LIBS = -lm -lpthread

TARGET = LWFQZip2
TARGET1 = LWMapping
TARGET2 = FQZip


all:	$(TARGET) $(TARGET1) $(TARGET2) 

$(TARGET):	$(OBJS)
	$(CC) -o $(TARGET) $(OBJS) $(LIBS)

$(TARGET1):	$(OBJS1)
	$(CC) -o $(TARGET1) $(OBJS1) $(LIBS)

$(TARGET2):	$(OBJS2)
	$(CXX) $(CXXFLAGS) $(OBJS2) -o $(TARGET2) $(LIBS)



$(OBJS): %.o: %.c
	$(CC) $(CFLAGS) -c $<

$(OBJS1): %.o: %.c
	$(CC) $(CFLAGS) -c $<


clean:
	-rm -f $(OBJS)
	-rm -f $(TARGET)
	-rm -f $(OBJS1)
	-rm -f $(TARGET1)
	-rm -f $(TARGET2)
