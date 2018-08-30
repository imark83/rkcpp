#ifndef __RK_COEFS_HPP__
#define __RK_COEFS_HPP__ value
// RK COEFICIENTS
#define A21	( 2.00000000000000E-01)
#define A31	( 7.50000000000000E-02)
#define A32	( 2.25000000000000E-01)
#define A41	( 9.77777777777778E-01)
#define A42	(-3.73333333333333E+00)
#define A43	( 3.55555555555556E+00)
#define A51	( 2.95259868922420E+00)
#define A52	(-1.15957933241884E+01)
#define A53	( 9.82289285169944E+00)
#define A54	(-2.90809327846365E-01)
#define A61	( 2.84627525252525E+00)
#define A62	(-1.07575757575758E+01)
#define A63	( 8.90642271774347E+00)
#define A64	( 2.78409090909091E-01)
#define A65	(-2.73531303602058E-01)
#define A71	( 9.11458333333333E-02)
#define A72	( 0.00000000000000E+00)		// zero!
#define A73	( 4.49236298292902E-01)
#define A74	( 6.51041666666667E-01)
#define A75	(-3.22376179245283E-01)
#define A76	( 1.30952380952381E-01)

#define C2	( 2.00000000000000E-01)
#define C3	( 3.00000000000000E-01)
#define C4	( 8.00000000000000E-01)
#define C5	( 8.88888888888889E-01)
#define C6	( 1.00000000000000E+00)
#define C7	( 1.00000000000000E+00)


#define B1	( 9.11458333333333E-02)
#define B2	( 0.00000000000000E+00)		// zero!
#define B3	( 4.49236298292902E-01)
#define B4	( 6.51041666666667E-01)
#define B5	(-3.22376179245283E-01)
#define B6	( 1.30952380952381E-01)
#define B7	( 0.00000000000000E+00)		// zero!


// COEFICIENTS OF ERROR ESTIMATION
#define E1	( 1.23263888888889E-03)
#define E2	( 0.00000000000000E+00)		// zero!
#define E3	(-4.25277029050614E-03)
#define E4	( 3.69791666666667E-02)
#define E5	(-5.08637971698113E-02)
#define E6	( 4.19047619047619E-02)
#define E7	(-2.50000000000000E-02)


#define CB1(th)	((((-1.00954861111111E+00*(th) + 2.88611111111111E+00)*(th) - 2.78541666666667E+00)*(th) + 1.00000000000000E+00)*(th))
#define CB2(th)	(0.00000000000000E+00)		// zero!
#define CB3(th)	(((2.27014076070680E+00*(th) - 5.60886492961965E+00)*(th) + 3.78796046720575E+00)*(th)*(th))
#define CB4(th) (((-2.16145833333333E+00*(th) + 4.50000000000000E+00)*(th) - 1.68750000000000E+00)*(th)*(th))
#define CB5(th) (((-1.32532429245283E+00*(th) + 1.26084905660377E+00)*(th) - 2.57900943396226E-01)*(th)*(th))
#define CB6(th)	(((2.22619047619048E+00*(th) - 3.03809523809524E+00)*(th) + 9.42857142857143E-01)*(th)*(th))

inline double denseEval(int nvar, double rkStage[], double x[], double step,
            int j, double th) {
  return x[j] + step*(CB1(th)*rkStage[j] + CB3(th)*rkStage[2*nvar+j] + CB4(th)*rkStage[3*nvar+j] + CB5(th)*rkStage[4*nvar+j] + CB6(th)*rkStage[5*nvar+j]);

}

#endif
