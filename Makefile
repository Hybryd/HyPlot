############
# MAKEFILE #
############

EXEC= hyplot
CC= g++
ERR_FLAGS= -W -Wall -pedantic -g
G_FLAGS= -lm -lpthread -lX11
C_FLAGS= $(G_FLAGS) -c
LD_FLAGS= $(G_FLAGS) -o
SRC= *.cpp
OBJ= $(SRC:.cpp=.o)
DOC= documentation

############

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) $^ $(LD_FLAGS) $@

.PHONY: clean mrproper

clean : $(OBJ)
	rm -f $^

mrproper: clean
	rm -f $(EXEC)

%.o : %.cpp
	$(CC) $(C_FLAGS) $^

create : $(DOC)
	doxygen -g $<

update : $(DOC)
	doxygen $<

doc : create update

