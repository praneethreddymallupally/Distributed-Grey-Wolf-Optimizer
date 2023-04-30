import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;



class DGWOImpl implements IDGWO{
	
	private int totalNumberOfCandidateSolutions;
	private double[][] candidateSolutions;
	private int totalNumberofDimensions;
	private double[] Lower;
	private double[] Upper;
	private double[] delta;
	private double[] beta;
	private double[] alpha;
	private fitnessclass iff;
	private int counter;
	private int maxiter;
	private double[] a;
	private double[] BESTVAL;
	private String identity;
	
	public DGWOImpl(int totalNumberOfCandidateSolutions, int numberOfTasks,  fitnessclass iff, int maxiter, double[] Lower, double[] Upper, String identity) {
		this.totalNumberOfCandidateSolutions = totalNumberOfCandidateSolutions;
		this.iff = iff;
		this.maxiter = maxiter;
		totalNumberofDimensions = numberOfTasks;
		alpha = new double[totalNumberofDimensions];
		beta = new double[totalNumberofDimensions];
		delta = new double[totalNumberofDimensions];
		counter = 0;
		a = new double[totalNumberofDimensions];
		BESTVAL  = new double[maxiter];
		candidateSolutions = new double[totalNumberOfCandidateSolutions][totalNumberofDimensions];
		this.Lower = Lower;
		this.Upper = Upper;
		this.identity = identity;
	}

	/* (non-Javadoc)
	 * @see IDGWO#sortAndIndex(double[][])
	 */
	public double[][] sortAndIndex(double[][] XXX) {
		
		Comparator<double[]> comparator = new Comparator<double[]>() {
            public int compare(double[] row1, double[] row2) {
            	double value1 = iff.fitnessfunction(row1);
            	double value2 = iff.fitnessfunction(row2);
                return Double.compare(value1, value2);
            }
        };
        Arrays.sort(XXX,comparator);
        return XXX;
	}

	/* (non-Javadoc)
	 * @see IDGWO#init()
	 */
	public void init() {
		for(int i=0; i<totalNumberOfCandidateSolutions; i++) {
			for(int j=0; j<totalNumberofDimensions; j++) {
				candidateSolutions[i][j] = (double)(Lower[j] + (Upper[j]- Lower[j])*Math.random());	
			}
		}
		candidateSolutions = sortAndIndex(candidateSolutions);
		for(int i=0; i<totalNumberofDimensions; i++) {
			alpha[i] = candidateSolutions[0][i];
		}
		for(int i=0; i<totalNumberofDimensions; i++) {
			beta[i] = candidateSolutions[1][i];
		}
		for(int i=0; i<totalNumberofDimensions; i++) {
			delta[i] = candidateSolutions[2][i];
		}
		System.out.println("+++++++++++ candidateSolutions +++++++++++");
		for(int i=0;i<totalNumberOfCandidateSolutions;i++) {
			for(int j=0; j<totalNumberofDimensions; j++) {
				System.out.print(candidateSolutions[i][j]);
				System.out.print(" ");
			}
			System.out.println();
		}
		System.out.println("+++++++++++ ++++++++++++++++++ +++++++++++");
	}

	/* (non-Javadoc)
	 * @see IDGWO#simplebounds(double[][])
	 */
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

	/* (non-Javadoc)
	 * @see IDGWO#solution()
	 */
	public double[][] solution() {
		double previous, after;
		if (counter == 0) {
	        init();
	    }
		System.out.println("Post Initialization....");
	    counter++;
	    int iter = 1;
	    while (iter < maxiter) {
	    	double worked=0, unworked=0;
	    	System.out.println("Iteration: "+iter+"Identity: "+identity);
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
					previous = X1[i][j];
					X1 = simplebounds(X1);
					after = X1[i][j];
					if (previous-after == 0) {
						worked+=1;
					}
					else {
						unworked+=1;
					}
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
					previous = X2[i][j];
					X2 = simplebounds(X2);
					after = X2[i][j];
					if (previous-after == 0) {
						worked+=1;
					}
					else {
						unworked+=1;
					}
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
					previous = X3[i][j];
					X3 = simplebounds(X3);
					after = X3[i][j];
					if (previous-after == 0) {
						worked+=1;
					}
					else {
						unworked+=1;
					}
	                candidateSolutions[i][j] = (X1[i][j] + X2[i][j] + X3[i][j]) / 3.0;

	            }
	        }
	        candidateSolutions = simplebounds(candidateSolutions);
	        candidateSolutions = sortAndIndex(candidateSolutions);

//	        for (int i = 0; i < totalNumberofDimensions; i++) {
//	        	candidateSolutions[totalNumberOfCandidateSolutions - 1][i] = candidateSolutions[0][i];
//	        }

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
//	        System.out.println("Worked = "+worked);
//	        System.out.println("Unworked = "+unworked);
	    }

	    double[][] out = new double[2][totalNumberofDimensions];
	    for (int i = 0; i < totalNumberofDimensions; i++) {
	        out[1][i] = alpha[i];
	    }
	    out[0][0] = iff.fitnessfunction(alpha);
	    return out;
	}

	/* (non-Javadoc)
	 * @see IDGWO#output()
	 */
	public double output() {
	    double[][] out = solution();
	    System.out.println("***********************************");
	    System.out.println("Optimized value = " + out[0][0]);
	      for(int i=0;i<totalNumberofDimensions;i++){
	    	  System.out.println("x["+i+"] = "+out[1][i]);
	      }
	    return out[0][0];
	}

	/* (non-Javadoc)
	 * @see IDGWO#migrateBest(int)
	 */
	public double[][] migrateBest(int n) {
		double migrants[][] = new double[n][totalNumberofDimensions];
	    for (int i = 0; i < n; i++) {
	        for (int j = 0; j < totalNumberofDimensions; j++) {
	            migrants[i][j] = candidateSolutions[i][j];
	        }
	    }
	    return migrants;
	}

	/* (non-Javadoc)
	 * @see IDGWO#replaceWorst(double[][])
	 */
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
	public double[] LOV(double[] arr) {
		double[] arrClone = Arrays.copyOf(arr, arr.length);
		
		HashMap < Double, Integer > map = new HashMap <Double, Integer > ();
		for (int j = 0; j < arr.length; j++) {
		    map.put(arr[j], j);
		}

		Arrays.sort(arrClone);
		
		double[] phi = new double[arr.length];
		for (int j = 0; j < arrClone.length; j++) {
			phi[map.get(arrClone[j])] = j;
		}

		double[] pi = new double[arr.length];
		for(int j=0; j< arrClone.length; j++) {
			pi[(int) phi[j]] = j;
		}
		return pi;
	}
}