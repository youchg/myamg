#################################################
HOME = ..
include $(HOME)/make.inc
#################################################
OBJSOURCE = $(wildcard *.c)
OBJLIB    = $(OBJSOURCE:%.c=%.o)

INC = -I$(CINC) $(INCEXTERN)
LIB = $(LIBAMG)
#################################################
.PHONY:info clean

lib: $(LIB) info

info:
	@echo " "
	@echo "Summary:"
	@echo "Libraries   are located in lib/"
	@echo " "

clean:
	@clear; clear;
	@$(RM) $(RMFLAG) *.o *~
	@$(RM) $(RMFLAG) $(LIBAMG)

#################################################
#  $@  目标文件
#  $^  所有的依赖文件
#  $<  第一个依赖文件
#################################################
%.o: %.c
	@echo " "
	@echo '$(CC) -c' $<
	@$(CC) $(CFLAG) -c $< $(INC)

$(LIB): $(OBJLIB)
	@echo " "
	@echo 'Creating library $(LIB)'
	@$(AR) $(ARFLAG) $(LIB) $^
#	@$(RANLIB) $(LIB)
#################################################
