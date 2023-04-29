import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;

interface DGWO{
	double[] getminval_index(double[] a);
	double[] getmaxval_index(double a[]);
	double[][] sort_and_index(double[][] XXX);
	void init();
	double[][] simplebounds(double s[][]);
	double[][] solution();
	double output();
	double[][] migrateBest(int n);
    void replaceWorst(double[][] migrants);
    void toStringnew();
}

abstract class f_xj
{
	abstract double fitnessfunction(double x[]);
}


class DGWOImpl implements DGWO{
	
	private int totalNumberOfCandidateSolutions;
	private double[][] candidateSolutions;
	private int totalNumberofDimensions;
	private int[] Lower;
	private int[] Upper;
	private double[] delta;
	private double[] beta;
	private double[] alpha;
	private f_xj iff;
	private int counter;
	private int maxiter;
	private double[] a;
	private double[] BESTVAL;
	
	public DGWOImpl() {
		
	}

	public double[] getminval_index(double[] a) {
		double m = 0.0;
	    double b[] = new double[a.length];
	    for (int i = 0; i < a.length; i++) {
	        b[i] = a[i];
	    }
	    double minval = a[0];
	    for (int i = 0; i < a.length; i++) {
	        if (a[i] < minval) {
	            minval = a[i];
	        }
	    }
	    for (int i = 0; i < a.length; i++) {
	        if (b[i] == minval) {
	            m = i;
	            break;
	        }
	    };
	    double[] dep = new double[2];
	    dep[0] = minval;
	    dep[1] = m;
	    return dep;
	}

	public double[] getmaxval_index(double[] a) {
	    double m = 0.0;
	    double b[] = new double[a.length];
	    for (int i = 0; i < a.length; i++) {
	        b[i] = a[i];
	    }
	    double maxval = a[0];
	    for (int j = 0; j < a.length; j++) {
	        if (a[j] > maxval) {
	            maxval = a[j];
	        }
	    }
	    for (int i = 0; i < b.length; i++) {
	        if (b[i] == maxval) {
	            m = i;
	            break;
	        }
	    }
	    double dep2[] = new double[2];
	    dep2[0] = maxval;
	    dep2[1] = m;
	    return dep2;
	}

	public double[][] sort_and_index(double[][] XXX) {
		double[] yval = new double[totalNumberOfCandidateSolutions];
	    for (int i = 0; i < totalNumberOfCandidateSolutions; i++) {
	        yval[i] = iff.fitnessfunction(XXX[i]);
	    }
	    ArrayList < Double > nfit = new ArrayList < Double > ();
	    for (int i = 0; i < totalNumberOfCandidateSolutions; i++) {
	        nfit.add(yval[i]);
	    }
	    ArrayList < Double > nstore = new ArrayList < Double > (nfit);
	    Collections.sort(nfit);
	    double[] ret = new double[nfit.size()];
	    Iterator < Double > iterator = nfit.iterator();
	    int ii = 0;
	    while (iterator.hasNext()) {
	        ret[ii] = iterator.next().doubleValue();
	        ii++;
	    }
	    int[] indexes = new int[nfit.size()];
	    for (int n = 0; n < nfit.size(); n++) {
	        indexes[n] = nstore.indexOf(nfit.get(n));
	    }
	    double[][] B = new double[totalNumberOfCandidateSolutions][totalNumberofDimensions];
	    for (int i = 0; i < totalNumberOfCandidateSolutions; i++) {
	        for (int j = 0; j < totalNumberofDimensions; j++) {
	            B[i][j] = XXX[indexes[i]][j];
	        }
	    }

	    return B;
	}

	public void init() {
		for(int i=0; i<totalNumberOfCandidateSolutions; i++) {
			for(int j=0; j<totalNumberofDimensions; j++) {
				candidateSolutions[i][j] = Lower[j] + (Upper[j]- Lower[j])*Math.random();
			}
		}
		candidateSolutions = sort_and_index(candidateSolutions);
		for(int i=0; i<totalNumberofDimensions; i++) {
			alpha[i] = candidateSolutions[0][i];
		}
		for(int i=0; i<totalNumberofDimensions; i++) {
			beta[i] = candidateSolutions[1][i];
		}
		for(int i=0; i<totalNumberofDimensions; i++) {
			delta[i] = candidateSolutions[2][i];
		}
	}

	public double[][] simplebounds(double[][] s) {
		for (int i = 0; i < totalNumberOfCandidateSolutions; i++) {
	        for (int j = 0; j < totalNumberofDimensions; j++) {
	            if (s[i][j] < Lower[j]) {
	                s[i][j] = Lower[j] * ((Upper[j] - Lower[j]) * Math.random());
	            }
	            if (s[i][j] > Upper[j]) {
	                s[i][j] = Lower[j] * ((Upper[j] - Lower[j]) * Math.random());
	            }
	        }
	    }
	    return s;
	}

	public double[][] solution() {
		if (counter == 0) {
	        init();
	    }
	    counter++;
	    int iter = 1;
	    while (iter < maxiter) {
	        for (int j = 0; j < totalNumberofDimensions; j++) {
	            a[j] = 2.0 - ((double) iter * (2.0 / (double) maxiter));
	        }

	        for (int i = 0; i < totalNumberOfCandidateSolutions; i++) {
	            for (int j = 0; j < totalNumberofDimensions; j++) {
	                // calculate X1
	            	double r1 = Math.random();
	                double r2 = Math.random();
	                double[] A1 = new double[totalNumberofDimensions];
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    A1[ii] = 2.0 * a[ii] * r1 - a[ii];
	                }
	                double[] C1 = new double[totalNumberofDimensions];
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    C1[ii] = 2.0 * r2;
	                }

	                double[][] X1 = new double[totalNumberOfCandidateSolutions][totalNumberofDimensions];
					X1[i][j] = alpha[j] - A1[j] * (Math.abs(C1[j] * alpha[j] - candidateSolutions[i][j]));
	                X1 = simplebounds(X1);
	                
	                // calculate X2
	                r1 = Math.random();
	                r2 = Math.random();
	                double[] A2 = new double[totalNumberofDimensions];
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    A2[ii] = 2.0 * a[ii] * r1 - a[ii];
	                }
	                double[] C2 = new double[totalNumberofDimensions];;
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    C2[ii] = 2.0 * r2;
	                }

	                double[][] X2 = new double[totalNumberOfCandidateSolutions][totalNumberofDimensions];;
					X2[i][j] = beta[j] - A2[j] * (Math.abs(C2[j] * beta[j] - candidateSolutions[i][j]));
	                X2 = simplebounds(X2);
	                
	                // calculate X3
	                r1 = Math.random();
	                r2 = Math.random();
	                double[] A3 = new double[totalNumberofDimensions];
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    A3[ii] = 2.0 * a[ii] * r1 - a[ii];
	                }
	                double[] C3 = new double[totalNumberofDimensions];;
					for (int ii = 0; ii < totalNumberofDimensions; ii++) {
	                    C3[ii] = 2.0 * r2;
	                }

	                double[][] X3 = new double[totalNumberOfCandidateSolutions][totalNumberofDimensions];;
					X3[i][j] = delta[j] - A3[j] * (Math.abs(C3[j] * delta[j] - candidateSolutions[i][j]));
	                X3 = simplebounds(X3);
	                
	                candidateSolutions[i][j] = (X1[i][j] + X2[i][j] + X3[i][j]) / 3.0;

	            }
	        }
	        candidateSolutions = simplebounds(candidateSolutions);
	        candidateSolutions = sort_and_index(candidateSolutions);

	        for (int i = 0; i < totalNumberofDimensions; i++) {
	        	candidateSolutions[totalNumberOfCandidateSolutions - 1][i] = candidateSolutions[0][i];
	        }

	        for (int i = 0; i < totalNumberofDimensions; i++) {
	            alpha[i] = candidateSolutions[0][i];
	        }
	        for (int i = 0; i < totalNumberofDimensions; i++) {
	            beta[i] = candidateSolutions[1][i];
	        }
	        for (int i = 0; i < totalNumberofDimensions; i++) {
	            delta[i] = candidateSolutions[2][i];
	        }

	        BESTVAL[iter] = iff.fitnessfunction(candidateSolutions[0]);

	        iter++;
	    }

	    double[][] out = new double[2][totalNumberofDimensions];
	    for (int i = 0; i < totalNumberofDimensions; i++) {
	        out[1][i] = alpha[i];
	    }
	    out[0][0] = iff.fitnessfunction(alpha);
	    return out;
	}

	public double output() {
	    double[][] out = solution();
	    System.out.println("Optimized value = " + out[0][0]);
	      for(int i=0;i<totalNumberofDimensions;i++){
	    	  System.out.println("x["+i+"] = "+out[1][i]);
	      }
	    return out[0][0];
	}

	public double[][] migrateBest(int n) {
		double migrants[][] = new double[n][totalNumberofDimensions];
	    for (int i = 0; i < n; i++) {
	        for (int j = 0; j < totalNumberofDimensions; j++) {
	            migrants[i][j] = candidateSolutions[i][j];
	        }
	    }
	    return migrants;
	}

	public void replaceWorst(double[][] migrants) {
		int n = candidateSolutions.length - 1;
	    int m = migrants.length - 1;
	    int c = 0;
	    while (c <= m) {
	        for (int i = n; i >= n - m; i--) {
	            for (int j = 0; j < totalNumberofDimensions; j++) {
	            	candidateSolutions[i][j] = migrants[c][j];
	            }
	        }
	        c++;
	    }
	}

	public void toStringnew() {
		double[][] in = solution();
	    System.out.println("Optimized value = " + in[0][0]);
	    for (int i = 0; i < totalNumberofDimensions; i++) {
	        System.out.println("x[" + i + "] = " + in[1][i]);
	    }
	}
	
}