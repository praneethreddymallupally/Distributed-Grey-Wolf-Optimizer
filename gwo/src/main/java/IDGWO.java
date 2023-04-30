
interface IDGWO {

	double[][] sortAndIndex(double[][] XXX);

	void init();

	double[][] simplebounds(double[][] s);

	double[][] solution();

	double output();

	double[][] migrateBest(int n);

	void replaceWorst(double[][] migrants);
	
	double[] LOV(double[] arr);
}