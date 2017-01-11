#################################################
HOME = ..
include $(HOME)/make.inc
#################################################
INC  = -I$(CINC)
LIB  = $(LIBAMG)
#################################################
SOURCE  = $(wildcard *.c)
TARGETS = $(SOURCE:%.c=%)
.PHONY: clean $(TARGETS) all
#################################################
test_list:test_list.exe
test_matrix:test_matrix.exe
test_cfsplit:test_cfsplit.exe
test_multigrid:test_multigrid.exe
test_cg:test_cg.exe
test_gs:test_gs.exe
test_linear_solver:test_linear_solver.exe
test_par_multigrid:test_par_multigrid.exe
test_par_cg:test_par_cg.exe
test_par_linear_solver:test_par_linear_solver.exe
test_par_linear_solver:test_par_linear_solver.exe

all: $(TARGETS)
#################################################
clean:
	@clear; clear;
	@$(RM) $(RMFLAG) *.o *~ *.exe
#################################################
%.exe: %.o
	@echo ""
	@echo '-------------------------------------------'
	@echo '$(CC) $(CFLAG) -o $@ $< $(LIB) $(LIBEXTERN)'
	@$(CC) $(CFLAG) -o  $@ $< $(LIB) $(LIBEXTERN)
	@echo '-------------------------------------------'
	@echo ""
%.o: %.c
	@echo '$(CC) $(CFLAG) -c $<'
	@$(CC) $(CFLAG) -c $< $(INC)
	@echo ""
#################################################
# wildcard 扩展通配符 -- 展开所有 .c 文件
#  $@  目标文件
#  $^  所有的依赖文件
#  $<  第一个依赖文件
#################################################
