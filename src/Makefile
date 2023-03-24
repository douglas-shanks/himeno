INCLUDES = -I. 

FC = ftn

FCFLAGS = -eF -homp -hacc

FCLINKFLAGS = -lcraymp

EXE_OMP=himeno_omp.exe

OBJS_OMP    = himeno_omp.o

EXE_ACC=himeno_acc.exe

OBJS_ACC    = himeno_acc.o

all: $(EXE_OMP) $(EXE_ACC)

$(EXE_OMP): $(OBJS_OMP)
	$(FC) $(FCLINKFLAGS) $< -o $@

$(EXE_ACC): $(OBJS_ACC)
	$(FC) $(FCLINKFLAGS) $< -o $@

%.o: %.f90
	$(FC) $(FCFLAGS) $(INCLUDES) -c $< -o $@	
	
himeno_omp.o: param.h

himeno_acc.o: param.h

clean:
	rm -rf $(EXE_OMP) $(OBJS_OMP) $(EXE_ACC) $(OBJS_ACC) *.i

cleanall: clean
	rm -rf my_output* *.rpt myrep* *.out *.lst *+pat *+apa *.ap2 *.xf *_exp.* *.acc.o *.acc.s


