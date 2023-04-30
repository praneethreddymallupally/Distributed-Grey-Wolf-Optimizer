
import java.io.File;
import java.util.Scanner;

class f1 extends fitnessclass //Gold stein f(x)=3.0 @x=(0,-1) -2<x[i]<2 i=1,2
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        double second = 0.0;
        first = (1.0 + (x[0] + x[1] + 1.0) * (x[0] + x[1] + 1.0) * (19.0 - 14.0 * x[0] + 3.0 * x[0] * x[0] - 14.0 * x[1] + 6.0 * x[0] * x[1] + 3.0 * x[1] * x[1]));
        second = 30.0 + (2.0 * x[0] - 3.0 * x[1]) * (2.0 * x[0] - 3.0 * x[1]) * (18.0 - 32.0 * x[0] + 12.0 * x[0] * x[0] + 48.0 * x[1] - 36.0 * x[0] * x[1] + 27 * x[1] * x[1]);
        return first * second;
    }
}

class f2 extends fitnessclass // Beale f(x)=0    @x=(3,0.5)   -4.5<x[i]<4.5, i = 1, 2.
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = ((1.5 - x[0] + x[0] * x[1]) * (1.5 - x[0] + x[0] * x[1])) + ((2.25 - x[0] + x[0] * x[1] * x[1]) * (2.25 - x[0] + x[0] * x[1] * x[1])) + ((2.625 - x[0] + x[0] * x[1] * x[1] * x[1]) * (2.625 - x[0] + x[0] * x[1] * x[1] * x[1]));
        return first;
    }
}

class f3 extends fitnessclass // Bohachecsky 1 f(x)=0  @x=(0.0,0.0)   -5.0<x[i]<5.0, i = 1, 2.
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = x[0] * x[0] + 2.0 * x[1] * x[1] - 0.3 * (Math.cos(Math.PI * 3.0 * x[0])) - 0.4 * Math.cos(4.0 * Math.PI * x[1]) + 0.7;
        return first;
    }
}

class f4 extends fitnessclass // Bohachecsky 2 f(x)=0  @x=(0.0,0.0)   -5.0<x[i]<5.0, i = 1, 2.
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = x[0] * x[0] + 2.0 * x[1] * x[1] - (0.3 * (Math.cos(Math.PI * 3.0 * x[0])) * Math.cos(4.0 * Math.PI * x[1])) + 0.3;
        return first;
    }
}

class f5 extends fitnessclass // Bohachecsky 3 f(x)=0  @x=(0.0,0.0)   -5.0<x[i]<5.0, i = 1, 2.
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = x[0] * x[0] + 2.0 * x[1] * x[1] - (0.3 * (Math.cos(Math.PI * 3.0 * x[0] + Math.PI * 4.0 * x[1]))) + 0.3;
        return first;
    }
}

class f6 extends fitnessclass // Booth  f(x)=0  @x=(1.0,3.0)   -10.0<x[i]<10.0, i = 1, 2.
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = (x[0] + 2.0 * x[1] - 7.0) * (x[0] + 2.0 * x[1] - 7.0) + (2.0 * x[0] + x[1] - 5.0) * (2.0 * x[0] + x[1] - 5.0);
        return first;
    }
}

class f7 extends fitnessclass // Branin  f(x)=0.397887  @x=(-pi,12.275),(pi,2.275),(9.42478,2.475)   -5.0<=x[0]<=10.0, 0.0<=x[1]<=15.0
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = ((x[1] - (5.1 * x[0] * x[0] / (4.0 * Math.PI * Math.PI)) + (5.0 * x[0] / Math.PI) - 6.0) * (x[1] - (5.1 * x[0] * x[0] / (4.0 * Math.PI * Math.PI)) + (5.0 * x[0] / Math.PI) - 6.0)) + (10.0 * (1.0 - (1.0 / (8.0 * 3.1415))) * Math.cos(x[0])) + 10.0;
        return first;
    }
}

class f8 extends fitnessclass // Colville  f(x)=0.0  @x=(1,1,1,1)   -10.0<=x[i]<=10.0 i=0,1,2,3
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = (100.0 * (x[0] - x[1] * x[1]) * (x[0] - x[1] * x[1])) + ((1.0 - x[0]) * (1.0 - x[0])) + (90.0 * (x[3] - x[2] * x[2]) * (x[3] - x[2] * x[2])) + ((1.0 - x[2]) * (1.0 - x[2])) + (10.1 * ((x[1] - 1.0) * (x[1] - 1.0) + (x[3] - 1.0) * (x[3] - 1.0))) + (19.8 * (x[1] - 1.0) * (x[3] - 1.0));
        return first;
    }
}

class f9 extends fitnessclass // Easom  f(x)=-1.0  @x=(pi,pi)   -100.0<=x[i]<=100.0 i=0,1
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = -Math.cos(x[0]) * Math.cos(x[1]) * Math.exp(-(x[0] - Math.PI) * (x[0] - Math.PI) - (x[1] - Math.PI) * (x[1] - Math.PI));
        return first;
    }
}

class f10 extends fitnessclass // Himmelblau f(x)=0.0  @x=(3.0,2.0),(-2.8051,3.1313),(-3.7793,-3.2831),(3.5844,-1.8481)   -6.0<=x[i]<=6.0 i=0,1
{
    double fitnessfunction(double x[]) {
        double first = 0.0;
        first = (((x[0] * x[0] + x[1] - 11.0) * (x[0] * x[0] + x[1] - 11.0)) + (x[0] + x[1] * x[1] - 7.0) * (x[0] + x[1] * x[1] - 7.0));
        return first;
    }
}

class f11 extends fitnessclass // Griewank f(x)=0.0  @x=(0,0)<---global minima     several local minimas      -600<x[i]<600 i=1,2,.. x.length 
{
    double fitnessfunction(double x[]) {
        double s = 0.0;
        double fact = 1.0;
        int m = x.length;
        for (int i = 0; i < m; i++) {
            s += x[i] * x[i];
        }
        for (int i = 0; i < m; i++) {
            fact *= Math.cos(x[i] / Math.sqrt(i + 1));
        }
        return (s / 4000.0) + 1.0 + (-fact);
    }
}

class f23 extends fitnessclass //Whitley's function         f(x)=0      @x=(0,0,0...) -10<x[i]<10 
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double s1 = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                s1 += (Math.pow((100.0 * (x[i] * x[i] - x[j]) * (x[i] * x[i] - x[j]) + (1.0 - x[j]) * (1.0 - x[j])), 2.0) / 4000.0) - Math.cos((100.0 * (x[i] * x[i] - x[j]) * (x[i] * x[i] - x[j]) + (1.0 - x[j]) * (1.0 - x[j]))) + 1.0;
            }
        }
        return s1;
    }
}

class f36 extends fitnessclass //Ackley�s function 2.9        f(x)=0;      @x=(0,0,0...)     -32.768<x[i]<32.768
{
    public double fitnessfunction(double x[]) {
        //��z�m� istenen fonksiyon	
        double a = 20.0;
        double b = 0.2;
        double c = 2. * Math.PI;

        int n = x.length;

        double r = Math.PI / 180;
        double top1 = 0.0;
        for (int i = 0; i < n; i++) {
            top1 += x[i] * x[i];
        }
        top1 = Math.sqrt(top1 / n);
        double top2 = 0.0;
        for (int i = 0; i < n; i++) {
            top2 += Math.cos(r * c * x[i]);
        }
        top2 = top2 / n;
        double top = -a * Math.exp(-b * top1) - Math.exp(top2) + a + Math.exp(1);
        return top;
    }
}

class f30 extends fitnessclass // Rosenbrock's valley     f(x)=0.0     -2.048<x[i]<2.048
{
    double fitnessfunction(double x[]) {
        //��z�m� istenen fonksiyon	
        int n = x.length;
        double ff = 0.0;
        for (int i = 0; i < n - 1; i++) {
            ff += (100.0 * (x[i + 1] - x[i] * x[i]) * (x[i + 1] - x[i] * x[i]) + (1.0 - x[i]) * (1.0 - x[i]));
        }
        return ff;
    }
}

class f63 extends fitnessclass // Generalized Schwefel�s Problem 2.26           
{
    public double fitnessfunction(double x[]) {
        int n = x.length;
        double top = 0.0;
        for (int i = 0; i < n; i++) {
            top += (x[i] * Math.sin(Math.sqrt(Math.abs(x[i]))));
        }
        return -top;
    }
}

class f101 extends fitnessclass // Alpine  -10.0<x[i]<10.0 f(x)=0  x(0,0,...,0)
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double t = 0.0;
        for (int i = 0; i < n; i++) {
            t += Math.abs(x[i] * Math.sin(x[i]) + 0.1 * x[i]);
        }
        return t;
    }
}

class f29 extends fitnessclass //Schaffer function         f(x)=0    @x=(0,0,...)  -100<x[i]<100  
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double s1 = 0.0;
        for (int i = 0; i < n; i++) {
            s1 += x[i] * x[i];
        }
        double s2 = Math.sqrt(s1);
        return 0.5 + ((Math.pow(Math.sin(s1), 2.0) - 0.5) / (1.0 + 0.001 * s1));
    }
}

class f34 extends fitnessclass //Rastrigin�s function 2.5        f(x)=0  @x=(0,0,...)     -5.12<x[i]<5.12
{
    public double fitnessfunction(double x[]) {
        //��z�m� istenen fonksiyon	
        double ff = 0;
        int n = x.length;
        for (int i = 0; i < n; i++) {
            ff += 10.0 + x[i] * x[i] - 10 * Math.cos(2.0 * Math.PI * x[i]);
        }
        return ff;
    }
}

class f1050 extends fitnessclass // dropwave  -5.12<x<5.12
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double sum = 0.0;
        for (int i = 0; i < n; i++) {
            sum += x[i] * x[i];
        }
        double nom = 1.0 + Math.cos(12.0 * Math.sqrt(sum));
        double sum1 = 0.0;
        for (int i = 0; i < n; i++) {
            sum1 += x[i] * x[i];
        }
        double denom = 0.5 * sum1 + 2.0;
        return -nom / denom;

    }
}

class f106 extends fitnessclass // Inverted Cosine Wave function,
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double t = 0.0;
        for (int i = 0; i < n - 1; i++) {
            t += Math.exp(-(x[i] * x[i] + x[i + 1] * x[i + 1] + 0.5 * x[i] * x[i + 1]) / 8.0) * Math.cos(4.0 * Math.sqrt(x[i] * x[i] + x[i + 1] * x[i + 1] + 0.5 * x[i] * x[i + 1]));
        }
        return -t;
    }
}

class f1021 extends fitnessclass //Levy function         f(x)=0   @x=(1,1,1...) -10<x[i]<10 
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double z[] = new double[n];
        for (int i = 0; i < n; i++) {
            z[i] = 1.0 + ((x[i] - 1.0) / 4.0);
        }
        double s = Math.pow(Math.sin(3.1415 * z[0]), 2.0);
        for (int i = 0; i < n - 1; i++) {
            s += Math.pow((z[i] - 1.0), 2.0) * (1.0 + 10.0 * Math.pow(Math.sin(3.1415 * z[i] + 1.0), 2.0));
        }
        return s + Math.pow(z[n - 1] - 1.0, 2.0) * (Math.pow(Math.sin(2.0 * 3.1415 * z[n - 1]), 2.0) + 1.0);
    }
}

class f1027 extends fitnessclass //Schwefel 2.22         f(x)=0    @x=(0,0,...)  -100<x[i]<100  
{
    double fitnessfunction(double x[]) {
        int n = x.length;
        double s1 = 0.0;
        double f1 = 1.0;
        for (int i = 0; i < n; i++) {
            s1 += Math.abs(x[i]);
            f1 *= Math.abs(x[i]);
        }
        return s1 + f1;
    }
}

class f1033 extends fitnessclass //Rotated hyper-ellipsoid function  -65.536<x[i]<65.536  f(x)=0   @x=(0,0,...)
{
    public double fitnessfunction(double x[]) {
        //��z�m� istenen fonksiyon	
        double ff = 0;
        int n = x.length;

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                ff += x[j] * x[j];
            }
        }
        return ff;
    }
}
class f1070 extends fitnessclass // shifted shphere function     f(x)=0  @x=(0,0,...)     -100<x[i]<100
{
    public double fitnessfunction(double x[]) {
        double[] shift = new double[x.length];
        double F_bias = -450, z = 0, lower = -100, upper = 100;
        double ff = 0;
        for (int j = 0; j < shift.length; j++) {
            shift[j] = lower + ((upper - lower) * Math.random());;
        }
        int n = x.length;
        for (int i = 0; i < n; i++) {
            z = x[i] - shift[i];
            ff += z * z;
        }
        return ff + F_bias;
    }
}
class f1071 extends fitnessclass // Shifted Schwefel’s problem2.2     f(x)=0  @x=(0,0,...)     -100<x[i]<100
{
    public double fitnessfunction(double x[]) {
        double[] shift = new double[x.length];
        double F_bias = -450, z = 0, lower = -100, upper = 100;
        double ff = 0;
        for (int j = 0; j < shift.length; j++) {
            shift[j] = lower + ((upper - lower) * Math.random());;
        }
        double sum0 = 0, sum1 = 0;
        int n = x.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                z = x[i] - shift[i];
                sum0 += z;
            }
            sum0 = sum0 * sum0;
            sum1 += sum0;
            sum0 = 0;
            z = 0;
        }
        return sum1 + F_bias;
    }
}

class f1073 extends fitnessclass // expanded function     f(x)=0  @x=(0,0,...)     -5<x[i]<5
{
    public double fitnessfunction(double x[]) {
        double[] sphere1 = new double[x.length];
        double F_bias = -130, lower = -5, upper = 5;
        double[] m_z = new double[x.length];;
        double[] m_o = new double[x.length];;
        for (int i = 0; i < m_o.length; i++) {
            m_o[i] -= lower + ((upper - lower) * Math.random());
        }
        for (int i = 0; i < m_o.length; i++) {
            m_z[i] = x[i] - m_o[i] + 1;
        }
        double result = 0;
        int n = x.length;

        for (int i = 1; i < x.length; i++) {
            result += F8(F2(m_z[i - 1], m_z[i]));
        }
        result += F8(F2(x[x.length - 1], x[0]));

        result = result + F_bias;

        return result;
    }
    static public double F8(double x) {
        return (((x * x) / 4000.0) - Math.cos(x) + 1.0);
    }
    static public double F2(double x, double y) {
        double temp1 = (x * x) - y;
        double temp2 = x - 1.0;
        return ((100.0 * temp1 * temp1) + (temp2 * temp2));
    }
}
class f1072 extends fitnessclass // shifted Rastrigin function     f(x)=0  @x=(0,0,...)     -5<x[i]<5
{
    public double fitnessfunction(double x[]) {
        double[] sphere1 = new double[x.length];
        double F_bias = -330, z = 0, lower = -5, upper = 5;
        double ff = 0;

        double result = 0;
        int n = x.length;
        for (int j = 0; j < sphere1.length; j++) {
            sphere1[j] = lower + ((upper - lower) * Math.random());;
        }

        for (int i = 0; i < n; i++) {
            z = x[i] - sphere1[i];
            result = result + ((Math.pow(z, 2)) - (10 * Math.cos(2 * Math.PI * z)) + 10);
        }

        result = result + F_bias;

        return result;
    }
}
public class DGWOTest {
    public static void main(String args[]) throws Exception {

        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f23 Whitley
        // double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0};

        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.120,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,};//  f34 Rastrigin
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12};

        // double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f29 schaffer
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,};

        //double[] Lower={-2.0,-2.0};// f1 gold stein
        //double[] Upper={2.0,2.0};

        //double[] Lower={-4.5,-4.5};// f2 beale
        //double[] Upper={4.5,4.5};

        //double[] Lower={-5.0,-5.0};// f3 bohac
        //double[] Upper={5.0,5.0};

        //double[] Lower={-5.0,-5.0};// f4 bohac
        //double[] Upper={5.0,5.0};

        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0};// f5 bohac
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0};

        // double[] Lower={-10.0,-10.0};// f6 booth
        // double[] Upper={10.0,10.0};

        // double[] Lower={-10.0,-10.0,-10.0,-10.0};// f8 collvile
        // double[] Upper={10.0,10.0,10.0,10.0};

        //double[] Lower={-6.0,-6.0};// f10 bohac
        //double[] Upper={6.0,6.0};

        //double[] Lower={-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0,-600.0};//  f11 griewank
        //double[] Upper={600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0,600.0};

        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f34 Rastrigin
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12};

        //double[] Lower={-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768};//  f36 ackley
        //double[] Upper={32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768};

        //double[] Lower={2.6,0.7,17.0,7.3,7.3,2.9,5.0}; //speed reducer
        //double[] Upper={3.6,0.8,28.0,8.3,8.3,3.9,5.5}; //f77

        //double[] Lower={-10.0,-10.0};//f699 
        //double[] Upper={10.0,10.0};//
        //double[] Lower={0.1,1.0,1.0,0.1}; //f76  1.7438 
        //double[] Upper={2.0,10.0,10.0,2.0};

        //double[] Lower={0.0,0.0};//f73 
        //double[] Upper={1.0,1.0};//

        //double[] Lower={0.0,0.25,2.00};//f75 
        //double[] Upper={2.0,1.30,15.00};//

        //double[] Lower={-2.0,-2.0,-2.0,-2.0,-2.0};   //const_f496--
        //double[] Upper={2.0,2.0,2.0,2.0,2.0};   

        //double[] Lower={5.49e-6,2.196e-3};   //const_f497--
        //double[] Upper={4.553,18.21 }; 

        //double[] Lower={0.0,0.0 };   //const_f498++
        //double[] Upper={10.0,10.0 }; 

        // double[] Lower={0.0,0.0 };   //const_f499++
        // double[] Upper={10.0,10.0 }; 

        //double[] Lower={-3.5,-3.5 };   //const_f500++
        //double[] Upper={2.5,2.5 }; 

        //double[] Lower={-5.0,-1.0,-5.0 };   //const_f501++
        //double[] Upper={5.0,3.0,5.0 }; 

        //double[] Lower={-2.0,1.0 };   //const_f502++
        //double[] Upper={2.0,6.0 }; 

        //double[] Lower={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0 };   //const_f503++
        //double[] Upper={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0 };

        //double[] Lower={0.0,0.0,0.0,0.0};   //const_f505++
        //double[] Upper={1.0,1.0,1.0,1.0};

        //double[] Lower={0.0,0.0,0.0};   //const_f506
        //double[] Upper={50.0,50.0,50.0};

        //double[] Lower={3.0,2.0,0.5};   //const_f507++
        //double[] Upper={5.0,4.0,2.0};

        //double[] Lower={-5.0,-5.0};   //const_f508++
        //double[] Upper={5.0,5.0};

        //double[] Lower={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};   //const_f509++
        //double[] Upper={1.0,1.0,1.0,1.0,1.0,1.0};

        //double[] Lower={0.25,1.5,};   //const_f510++
        //double[] Upper={1.0,2.0*Math.PI};

        //double[] Lower={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};   //const_f511++
        //double[] Upper={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

        //double[] Lower={-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};   //const_f512++
        //double[] Upper={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};

        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0};   //const_f513--
        //double[] Upper={5.0,5.0,5.0,5.0,5.0}; 

        //double[] Lower={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,-3.0};   //const_f514--
        //double[] Upper={3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0,3.0}; 

        //double[] Lower={0.0,0.0};   //const_f515++
        //double[] Upper={1.0,1.0};

        //double[] Lower={-3.0,-3.0,-3.0,-3.0,-3.0,-3.0,};  //const_f516-- 
        //double[] Upper={3.0,3.0,3.0,3.0,3.0,3.0,};

        //double[] Lower={0.0,0.0};   //const_f517++
        //double[] Upper={2.0*Math.PI,2.0*Math.PI};

        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0};   //const_f518--
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0}; 

        //double[] Lower={-5.0,-5.0};   //const_f519
        //double[] Upper={5.0,5.0};

        //double[] Lower={0.06,0.06,0.06};   //const_f520++
        //double[] Upper={1.0,1.0,1.0};

        //double[] Lower={-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048};//  f30 Rosenbrock
        //double[] Upper={2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048};

        //Unimodal
        // f31 Sphere
        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f31 sphere
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12 };
        // f27 Schfewel 2.22
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f11 griewank
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,};
        // f30 Rosenbrock //
        //double[] Lower={-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048};//  f30 Rosenbrock
        //double[] Upper={2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048};
        //f40 step 
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f40 step
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        // f25 quartic//
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0};//  f11 griewank
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0 };
        // f63 Generalized Schwefel�s Problem 2.26
        //double[] Lower={-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0,-500.0};//  f11 griewank
        //double[] Upper={500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0,500.0 };
        // f34 Rastrigin
        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f34 Rastrigin
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12};
        //36 ackley
        //double[] Lower={-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768};//  f36 ackley
        //double[] Upper={32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768};
        // f11 griewank
        // double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f11 griewank
        // double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0};
        //f41 penal1
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f41 Penalized
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        //f42 penal2
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f41 Penalized
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        // f101 alpine//
        // double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0};//  f11 griewank
        // double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0 };
        //f106 inverted cosine
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f40 step
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        //f107 pathologic
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f40 step
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        //f108 Salomon
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f40 step
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
        //f19 zakharov
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0};//  f19 Zakharov
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0};                  
        //f20 levy                
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f20 levy
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,};

        //////////// 50D
        //f20 levy                
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,};//  f20 levy
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,};
        //f40 step 
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f40 step
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,};
        //f41 penal1
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f41 Penalized
        //double[] Upper={5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0,};
        //f19 zakharov//
        //double[] Lower={-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,};//  f19 Zakharov
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,};                  
        //f36 ackley //
        //double[] Lower={-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768,-32.768};//  f36 ackley
        //double[] Upper={32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768,32.768};
        //f8 collvile
        //double[] Lower={-10.0,-10.0,-10.0,-10.0}; //f8 colville
        //double[] Upper={10.0,10.0,10.0,10.0};   
        // f11 griewank
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0};//  f11 griewank
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
        // f34 Rastrigin //
        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f34 Rastrigin
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12};
        // f30 Rosenbrock //
        //double[] Lower={-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.0480,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048,-2.048};//  f30 Rosenbrock
        //double[] Upper={2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048,2.048, 2.048, 2.048, 2.0480, 2.048, 2.048, 2.048, 2.048, 2.048, 2.048, 2.048};
        // f31 Sphere//
        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f31 sphere
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12 };
        // f37 Dropwave//
        //double[] Lower={-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12,-5.12};//  f31 sphere
        //double[] Upper={5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12,5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12, 5.12 };
        // f80 Schwefel
        //double[] Lower={-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500,-500};//  f31 sphere
        //double[] Upper={500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500,500, 500, 500, 500, 500, 5000, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500 };
        // f25 quartic//
        //double[] Lower={-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,-10.0};//  f11 griewank
        //double[] Upper={10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};

        int size = 10; // number of decision variables of the problem
        double[] Lower = {0,1,2,3,4,5,4,3,2,1};
        double[] Upper = {6,7,8,9,9,9,9,8,7,6};
        f101 ff = new f101();
        int Maxiter = 1000;
        int Fm = 50; // migration frequency
        int N = 60; // population size
        int numSubpopulation = 10; // number of islands
        double Rm = 0.3; // migration 	rate
        int numberofMigrants = (int)(N / numSubpopulation * Rm); //number of migrants

        IDGWO wo0 = new DGWOImpl(N, size, ff, Maxiter, Lower, Upper,"total");
        
        IDGWO wo1 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-0");
        IDGWO wo2 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-1");
        IDGWO wo3 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-2");
        IDGWO wo4 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-3");
        IDGWO wo5 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-4");
        IDGWO wo6 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-5");
        IDGWO wo7 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-6");
        IDGWO wo8 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-7");
        IDGWO wo9 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-8");
        IDGWO wo10 = new DGWOImpl(N / numSubpopulation, size, ff, Fm, Lower, Upper,"Island-9");

        wo0.output();
        double bestOptimizedValue = 100000000, result1, result2, result3, result4, result5, result6, result7, result8, result9, result10;

        for (int i = 0; i < Maxiter / Fm; i++) {
        	System.out.println("No. of Migrations completed until now: "+i);
            result1 = wo1.output();
            if (result1 < bestOptimizedValue)
                bestOptimizedValue = result1;
            	
            result2 = wo2.output();
            if (result2 < bestOptimizedValue)
                bestOptimizedValue = result2;
            result3 = wo3.output();
            if (result3 < bestOptimizedValue)
                bestOptimizedValue = result3;
            result4 = wo4.output();
            if (result4 < bestOptimizedValue)
                bestOptimizedValue = result4;
            result5 = wo5.output();
            if (result5 < bestOptimizedValue)
                bestOptimizedValue = result5;
            result6 = wo6.output();
            if (result6 < bestOptimizedValue)
                bestOptimizedValue = result6;
            result7 = wo7.output();
            if (result7 < bestOptimizedValue)
                bestOptimizedValue = result7;
            result8 = wo8.output();
            if (result8 < bestOptimizedValue)
                bestOptimizedValue = result8;
            result9 = wo9.output();
            if (result9 < bestOptimizedValue)
                bestOptimizedValue = result9;
            result10 = wo10.output();
            if (result10 < bestOptimizedValue)
                bestOptimizedValue = result10;

            if (Fm % 10 == 0 || Fm % 10 == 1) {
                wo2.replaceWorst(wo1.migrateBest(numberofMigrants));
                wo3.replaceWorst(wo2.migrateBest(numberofMigrants));
                wo4.replaceWorst(wo3.migrateBest(numberofMigrants));
                wo5.replaceWorst(wo4.migrateBest(numberofMigrants));
                wo6.replaceWorst(wo5.migrateBest(numberofMigrants));
                wo8.replaceWorst(wo7.migrateBest(numberofMigrants));
                wo9.replaceWorst(wo8.migrateBest(numberofMigrants));
                wo10.replaceWorst(wo9.migrateBest(numberofMigrants));
                wo7.replaceWorst(wo10.migrateBest(numberofMigrants));
                wo1.replaceWorst(wo7.migrateBest(numberofMigrants));
            } else if (Fm % 10 == 02 || Fm % 10 == 3) {
                wo3.replaceWorst(wo1.migrateBest(numberofMigrants));
                wo5.replaceWorst(wo3.migrateBest(numberofMigrants));
                wo4.replaceWorst(wo5.migrateBest(numberofMigrants));
                wo6.replaceWorst(wo4.migrateBest(numberofMigrants));
                wo8.replaceWorst(wo6.migrateBest(numberofMigrants));
                wo7.replaceWorst(wo8.migrateBest(numberofMigrants));
                wo10.replaceWorst(wo7.migrateBest(numberofMigrants));
                wo2.replaceWorst(wo10.migrateBest(numberofMigrants));
                wo9.replaceWorst(wo2.migrateBest(numberofMigrants));
                wo1.replaceWorst(wo9.migrateBest(numberofMigrants));
            } else if (Fm % 10 == 4 || Fm % 10 == 5) {
                wo1.replaceWorst(wo10.migrateBest(numberofMigrants));
                wo2.replaceWorst(wo1.migrateBest(numberofMigrants));
                wo3.replaceWorst(wo2.migrateBest(numberofMigrants));
                wo4.replaceWorst(wo3.migrateBest(numberofMigrants));
                wo6.replaceWorst(wo4.migrateBest(numberofMigrants));
                wo5.replaceWorst(wo6.migrateBest(numberofMigrants));
                wo8.replaceWorst(wo5.migrateBest(numberofMigrants));
                wo7.replaceWorst(wo8.migrateBest(numberofMigrants));
                wo9.replaceWorst(wo7.migrateBest(numberofMigrants));
                wo10.replaceWorst(wo9.migrateBest(numberofMigrants));
            } else {
                wo5.replaceWorst(wo6.migrateBest(numberofMigrants));
                wo1.replaceWorst(wo5.migrateBest(numberofMigrants));
                wo3.replaceWorst(wo1.migrateBest(numberofMigrants));
                wo7.replaceWorst(wo3.migrateBest(numberofMigrants));
                wo10.replaceWorst(wo7.migrateBest(numberofMigrants));
                wo2.replaceWorst(wo10.migrateBest(numberofMigrants));
                wo9.replaceWorst(wo2.migrateBest(numberofMigrants));
                wo8.replaceWorst(wo9.migrateBest(numberofMigrants));
                wo4.replaceWorst(wo8.migrateBest(numberofMigrants));
                wo6.replaceWorst(wo4.migrateBest(numberofMigrants));
            }

        }
        System.out.println("***************************");
        System.out.println("***************************");
        System.out.println("***************************");
        System.out.println("Optimized value12 = " + bestOptimizedValue);
        System.out.println("***************************");
        System.out.println("***************************");
        System.out.println("***************************");

    }

}