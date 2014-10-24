import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Scanner;

class absorb{

	public static void main(String[] arg) throws IOException{
	double nu=1., ltot=1., la2=1., lc2=1., sindex=1., nu_old=1.0, ltot_old=1.0;
	double zz = 0.6035, tau, exptau;
	FileReader fin = new FileReader("emission.dat");
	FileWriter fout = new FileWriter("EBLcorrected.dat");
	Scanner src = new Scanner(fin);
	while (src.hasNextDouble()) {
		nu = src.nextDouble();
		ltot = src.nextDouble();
		la2 = src.nextDouble();
		lc2 = src.nextDouble();
		sindex = src.nextDouble();

		tau = franceschini.tau_IRA(nu,zz);
		exptau = Math.exp(-tau);
		ltot *= exptau;
		la2 *= exptau;
		lc2 *= exptau;
		sindex = 0.;
		if(ltot>1e-20)sindex = Math.log(ltot/ltot_old)/Math.log(nu/nu_old);

		nu_old = nu;
		ltot_old = ltot;
		fout.write(String.format("%14.7e%14.7e%14.7e%14.7e%14.6e\n",
		Math.max(nu,1e-60),Math.max(ltot,1e-60),Math.max(la2,1e-60),Math.max(lc2,1e-60),sindex));

		}
	fin.close();
	fout.close();
	}
}
