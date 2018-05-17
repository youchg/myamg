#################################################
HOME = ..
include $(HOME)/../make.inc
#################################################
INC  = -I$(CINC) $(INCEXTERN)
LIB  = $(LIBAMG)
#################################################
SOURCE  = $(wildcard *.c)
TARGETS = $(SOURCE:%.c=%)
.PHONY: clean $(TARGETS)
#################################################
ex_linear_solver:ex_linear_solver.exe
ex_linear_solver_shift:ex_linear_solver_shift.exe

ex_eigen_solver :ex_eigen_solver.exe
ex_par_eigen_solver :ex_par_eigen_solver.exe

ex1_poisson     :ex1_poisson.exe
ex2_lshape      :ex2_lshape.exe
ex3_1000_1_001_1:ex3_1000_1_001_1.exe
ex4_10_1_10_1   :ex4_10_1_10_1.exe
ex5_4_1_4_1     :ex5_4_1_4_1.exe

ex6_poisson_delaunay:ex6_poisson_delaunay.exe
ex7_lshape_delaunay:ex7_lshape_delaunay.exe
ex8_lshape_delaunay_adapt:ex8_lshape_delaunay_adapt.exe

ex9_hydrogen      :ex9_hydrogen.exe
ex10_fem3dlaplace :ex10_fem3dlaplace.exe
ex11_fem3dhydrogen:ex11_fem3dhydrogen.exe

ex12_fem2d_anisotropic:ex12_fem2d_anisotropic.exe

ex13_jiusuo:ex13_jiusuo.exe

ex15_poisson_augmented:ex15_poisson_augmented.exe

ex16_poisson_augmented_one:ex16_poisson_augmented_one.exe
#################################################
clean:
	@clear; clear;
	@$(RM) $(RMFLAG) *.o *~ *.exe
#################################################
%.exe: %.o $(LIB)
	@echo '$(CC) $(CFLAG) -o $@ $< $(LIB) $(LIBEXTERN)'
	@$(CC) $(CFLAG) -o  $@ $< $(LIB) $(LIBEXTERN)
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
