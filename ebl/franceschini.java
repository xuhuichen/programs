public class franceschini{

/**************************************************************************
*  IR absorbtion 
*
*  franceschini 2008
*
*
****************************************************************************/
final static int ZDIM_Franceschini = 9;
final static int EDIM_Franceschini = 50;
final static int MDIM_Franceschini = 450;

static double[] Z3 = {0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0} ;

static double[] E3 = 
{0.0200,0.0240,0.0289,0.0347,0.0417,0.0502,0.0603,0.0726,0.0873,0.104 ,0.126 ,0.151, 
0.182,0.219 ,0.263 ,0.316 ,0.381 ,0.458 ,0.550 ,0.662 ,0.796 ,0.957 ,1.15  ,1.38 
,1.66  ,2.000 ,2.40  ,2.89  ,3.47  ,4.17  ,5.02  ,6.03  ,7.26  ,8.73  ,10.5  ,12.6 
,15.2  ,18.2  ,21.9  ,26.4  ,31.7  ,38.1  ,45.8  ,55.1  ,66.2  ,79.6  ,95.7  ,115. 
,138.  ,166.   } ; 

static double[] T3 = {0.0000, 0.0000, 0.0000, 0.000000, 0.0000021, 
0.004933, 0.0399, 0.1157, 0.2596, 
0.0000, 0.0000, 0.0000, 0.000000, 0.000188, 0.01284, 0.0718, 0.1783, 0.3635, 
0.0000, 0.0000, 0.0000, 0.000000, 0.001304, 0.0279, 0.1188, 0.2598, 0.4919, 
0.0000, 0.0000, 0.0000, 0.000488, 0.004558, 0.0533, 0.1833, 0.3635, 0.6517, 
0.0000, 0.0000, 5.254E-05, 0.002276, 0.01157, 0.0921, 0.2689, 0.4967, 0.8548, 
0.0000, 9.445E-05, 5.408E-04, 0.006575, 0.02436, 0.1480, 0.3836, 0.6745, 1.118, 
1.0976E-04, 4.241E-04, 1.915E-03, 0.014592, 0.04512, 0.2275, 0.5434, 0.9179, 1.465, 
3.0882E-04, 1.103E-03, 4.548E-03, 0.02771, 0.07684, 0.3430, 0.7707, 1.251, 1.917, 
6.5619E-04, 2.258E-03, 8.903E-03, 0.04808, 0.1248, 0.5137, 1.092, 1.703, 2.503, 
1.2130E-03, 4.097E-03, 1.582E-02, 0.07958, 0.1984, 0.7640, 1.537, 2.302, 3.249, 
2.1063E-03, 7.039E-03, 2.685E-02, 0.1284, 0.3109, 1.120, 2.133, 3.073, 4.181, 
3.5291E-03, 1.167E-02, 4.406E-02, 0.2031, 0.4780, 1.607, 2.905, 4.042, 5.318, 
5.7051E-03, 1.872E-02, 7.010E-02, 0.3134, 0.7163, 2.247, 3.875, 5.225, 6.673, 
8.9183E-03, 2.907E-02, 0.1082, 0.4696, 1.040, 3.056, 5.055, 6.627, 8.241, 
1.3517E-02, 4.378E-02, 0.1618, 0.6809, 1.461, 4.042, 6.438, 8.226, 9.997, 
1.9793E-02, 6.367E-02, 0.2338, 0.9517, 1.981, 5.192, 7.989, 9.977, 11.89, 
2.7938E-02, 8.935E-02, 0.3256, 1.281, 2.594, 6.474, 9.650, 11.81, 13.89, 
3.7957E-02, 0.1205, 0.4356, 1.661, 3.284, 7.836, 11.34, 13.67, 15.93, 
4.9558E-02, 0.1563, 0.5607, 2.082, 4.023, 9.214, 13.01, 15.51, 18.08, 
6.2291E-02, 0.1953, 0.6961, 2.524, 4.779, 10.55, 14.63, 17.39, 20.45, 
7.5753E-02, 0.2364, 0.8373, 2.967, 5.517, 11.82, 16.25, 19.49, 23.27, 
8.9194E-02, 0.2768, 0.9750, 3.389, 6.210, 13.03, 18.04, 22.02, 26.81, 
0.1019, 0.3152, 1.105, 3.779, 6.846, 14.29, 20.21, 25.22, 31.33, 
0.1136, 0.3501, 1.223, 4.129, 7.432, 15.73, 22.98, 29.37, 37.23, 
0.1240, 0.3810, 1.327, 4.444, 8.010, 17.54, 26.58, 34.78, 45.09, 
0.1329, 0.4076, 1.419, 4.747, 8.652, 19.87, 31.31, 41.95, 55.80, 
0.1409, 0.4318, 1.504, 5.079, 9.452, 22.96, 37.67, 51.72, 70.71, 
0.1486, 0.4560, 1.596, 5.498, 10.52, 27.08, 46.30, 65.17, 92.14, 
0.1579, 0.4863, 1.714, 6.075, 11.96, 32.66, 58.24, 84.48, 124.0, 
0.1710, 0.5284, 1.879, 6.875, 13.92, 40.39, 75.45, 113.2, 172.1, 
0.1896, 0.5887, 2.113, 7.952, 16.57, 51.39, 101.3, 157.7, 245.9, 
0.2162, 0.6732, 2.431, 9.421, 20.24, 67.70, 141.4, 227.3, 357.0, 
0.2512, 0.7847, 2.858, 11.45, 25.57, 92.73, 204.9, 335.3, 519.3, 
0.3017, 0.9447, 3.464, 14.41, 33.52, 132.2, 304.3, 496.4, 747.0, 
0.3732, 1.171, 4.334, 18.87, 45.81, 195.0, 454.1, 724.1, 1048., 
0.4795, 1.513, 5.663, 25.76, 65.21, 292.5, 666.6, 1027., 1426., 
0.6455, 2.048, 7.723, 36.60, 95.98, 435.4, 948.5, 1407., 1873., 
0.8984, 2.871, 10.93, 53.79, 143.7, 630.5, 1299., 1852., 2372., 
1.297, 4.162, 15.99, 80.41, 214.3, 878.2, 1705., 2339., 2897., 
1.917, 6.177, 23.86, 119.9, 311.5, 1172., 2145., 2843., 3412., 
2.856, 9.181, 35.47, 174.8, 435.3, 1494., 2587., 3323., 3887., 
4.211, 13.47, 51.78, 245.0, 580.9, 1824., 3004., 3752., 4303., 
6.038, 19.14, 72.76, 327.4, 739.9, 2134., 3362., 4102., 4618., 
8.285, 26.00, 97.51, 416.9, 899.0, 2403., 3638., 4352., 4835., 
10.82, 33.59, 124.3, 506.3, 1045., 2616., 3825., 4499., 4940., 
13.48, 41.42, 151.2, 587.7, 1169., 2756., 3915., 4540., 4945., 
16.04, 48.81, 175.8, 655.4, 1263., 2823., 3915., 4540., 4945., 
18.24, 54.98, 195.8, 705.7, 1320., 2823., 3915., 4540., 4945., 
20.01, 59.82, 210.8, 735.5, 1340., 2823., 3915., 4540., 4945., 
21.20, 63.03, 219.8, 744.0, 1340., 2823., 3915., 4540., 4945.};

static double Interpolate3(double ee, double zz, boolean print_flag) {
      int ei = 1, zi = 1;
      int col, row;
      
      double x1=0, y1=0, x2=0, y2=0;
      double z11=0, z12=0, z21=0, z22=0;
      double a, b, c, d, t;

      if ((ee >= E3[0]) && (ee <= E3[EDIM_Franceschini-1]) &&
          (zz >= Z3[0]) && (zz <= Z3[ZDIM_Franceschini-1])) {

        if (print_flag) System.out.println( "Input: E=%7.5f z=%7.5f\n\n" +ee+ "," +zz);
      
        for (col=0; col <= ZDIM_Franceschini-2; col++) {
           if ((zz >= Z3[col]) && (zz <= Z3[col+1])) {
    	     zi = col;
	     break;
	   }
        }
            
        for (row=0; row <= EDIM_Franceschini-2; row++) {
           if ((ee >= E3[row]) && (ee <= E3[row+1])) {
	     ei = row;
	     break;
	   }
        }
            
        if (print_flag) System.out.println( "Found: E_min=%7.5f E_max=%7.5f ei=%d\n" +E3[ei]+ "," +E3[ei+1]+ ","+ei);
        if (print_flag) System.out.println( "Found: z_min=%7.5f z_max=%7.5f zi=%d\n\n" +Z3[zi]+ "," +Z3[zi+1]+ "," +zi);
      
        x1  = Z3[zi];
        x2  = Z3[zi+1];
        y1  = E3[ei];
        y2  = E3[ei+1];
      
        z11 = T3[ZDIM_Franceschini*ei+zi];           // [ei][zi];
        z12 = T3[ZDIM_Franceschini*(ei+1)+zi];       // [ei+1][zi];
        z21 = T3[ZDIM_Franceschini*ei+zi+1];         // [ei][zi+1];
        z22 = T3[ZDIM_Franceschini*(ei+1)+zi+1];     // [ei+1][zi+1];
      
        a = (z11 - z12 - z21 + z22)             / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        b = (y1*z12 - y2*z11 - y1*z22 + y2*z21) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        c = (x1*z21 - x2*z11 - x1*z22 + x2*z12) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
        d = (x1*y1*z22 - x1*y2*z21 - x2*y1*z12 + x2*y2*z11) / (x1*y1 - x1*y2 - x2*y1 + x2*y2);
      
        // the function
        t = a*ee*zz + b*zz + c*ee + d;
      } else {
        t = 0;
      }
      
      if (print_flag) System.out.println( "(x1=%7.5f, y1=%7.5f)                  (x2=%7.5f, y1=%7.5f)\n" +x1+ "," +y1+ "," +x2+ "," +y1);
      if (print_flag) System.out.println( "                    (x=%7.5f, y=%7.5f)\n" +zz+ "," +ee);
      if (print_flag) System.out.println( "(x1=%7.5f, y2=%7.5f)                  (x2=%7.5f, y2=%7.5f)\n\n" +x1+ "," +y2+ "," +x2+ "," +y2);
     
    if (print_flag) System.out.println( "(z11=%7.3f, z11_dim=%d)           (z21=%7.3f, z12_dim=%d)\n" +z11+ "," +(ZDIM_Franceschini*ei+zi)+ "," +z21+ "," +(ZDIM_Franceschini*ei+zi+1) );
      if (print_flag) System.out.println( "             (z=%7.3f)\n" +t);
      if (print_flag) System.out.println( "(z12=%7.3f, z12_dim=%d)           (z22=%7.3f, z12_dim=%d)\n" +z12+ "," +(ZDIM_Franceschini*(ei+1)+zi)+ "," +z22+ "," +(ZDIM_Franceschini*(ei+1)+zi+1));
            
      return t;
}


static double tau_IRA(double nu, double zz) {
      double E_keV, E_GeV, E_TeV, tau;
      double keV = 2.417898940e17;
      
      E_keV = nu / keV;
      E_GeV = E_keV / 1.0e+6;
      E_TeV = E_keV / 1.0e+9;
      
      tau = Interpolate3(E_TeV, zz, false);     
      return tau;
}

static double[] tau_IRA(double[] nu, double zz){
	double[] tau = new double[nu.length];

	for(int i=0; i < nu.length; i++){
	tau[i]=tau_IRA(nu[i],zz);
	}
	return tau;
}
}

