LIB_OBJ = libinformation.o libbasicanalysis.o libadvancedanalysis.o librfifiltering.o libshm.o libacquiredata.o 


all: $(LIB_OBJ)
	ar rc libajax.a $^
	@echo "Done compiling libraries"
%.o : 
	$(CC) $(@:.o=.cpp) -o $@ $(CFLAGS) 
	
clean:
	rm -rf *.o
	rm -rf *.a
