# an example makefile for pthreads GOL
# 
#CC = gcc
CC = nvcc
# for a build with debug info in executable:
# CFLAGS = -g -Wall -Werror=vla 
CFLAGS = -g -G 
AR = /usr/bin/ar -rvs 


# for a build with compiler optimization:
#CFLAGS = -O2 -Wall -Werror=vla

LIBS = -lm -lpthread 

DIR=lib

INCLUDES=$(DIR)/eig3.h $(DIR)/eig3c.h

TARGET1=s_zipper
TARGET2=p_zipper
TARGET3=p_zipper_old

all: $(TARGET1) $(TARGET2) #$(TARGET3)

#all: $(TARGET2)

$(DIR)/myopengllib.o: $(DIR)/myopengllib.cu 
	${CC} ${CFLAGS} -c $(DIR)/myopengllib.cu -lglut -lGL -lGLU -lGLEW 

$(DIR)/libmyopengllib.a: $(DIR)/myopengllib.o
	${AR} $(DIR)/libmyopengllib.a  $(DIR)/myopengllib.o 

$(DIR)/eig3.o: $(DIR)/eig3.c $(INCLUDES)
	$(CC) $(CFLAGS) -c $(DIR)/eig3.c $(LIBS)

$(DIR)/eig3c.o: $(DIR)/eig3c.cpp $(INCLUDES)
	$(CC) $(CFLAGS) -c $(DIR)/eig3c.cpp $(LIBS)


$(TARGET1): $(DIR)/eig3.o $(TARGET1).c $(INCLUDES)
	$(CC) $(CFLAGS) -o $(TARGET1) $(DIR)/eig3.o $(TARGET1).c $(LIBS) 

$(TARGET2): $(DIR)/eig3c.o $(TARGET2).cu $(INCLUDES) $(DIR)/libmyopengllib.a
	$(CC) $(CFLAGS) -o $(TARGET2) $(DIR)/eig3c.o $(TARGET2).cu $(LIBS) -L./lib -lmyopengllib -lglut -lGLU -lGLEW


$(TARGET3): $(DIR)/eig3c.o $(TARGET3).cu $(INCLUDES) $(DIR)/libmyopengllib.a
	$(CC) $(CFLAGS) -o $(TARGET3) $(DIR)/eig3c.o $(TARGET3).cu $(LIBS) -L./lib -lmyopengllib -lglut -lGLU -lGLEW

clean:
	$(RM) $(TARGET1) $(TARGET2) $(DIR)/myopengllib.o $(DIR)/libmyopengllib.a	


