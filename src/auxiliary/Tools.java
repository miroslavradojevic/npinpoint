package auxiliary;

import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.*;
import java.io.File;
import java.io.FilenameFilter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

/** various computational tools and file management tools */
public class Tools {

    static double VERY_SMALL_POSITIVE = 0.000001;
    static float  VERY_SMALL_POSITIVE_FLOAT = 0.000001f;

    public static ImagePlus convertToFloatImage(
            ImagePlus inim
    )
    {

        int W = inim.getWidth();
        int H = inim.getHeight();

        ImageProcessor ip = new FloatProcessor(W, H);
        for (int i = 0; i < H*W; i++ ) {
            ip.setf(i, inim.getProcessor().getPixelValue(i%W, i/W));
        }

        ImagePlus outim = new ImagePlus(inim.getTitle()+"ToFloat", ip);
        return outim;

    }

    public static String getExtension(String filePath)
    {
        String extension = "";
        int dotIdx = filePath.lastIndexOf('.');
        if (dotIdx > 0)
            extension = filePath.substring(dotIdx+1);
        return  extension;
    }

    public static String removeExtension(String filePath)
    {
        String extension = "";
        int dotIdx = filePath.lastIndexOf('.');
        if (dotIdx > 0)
            extension = filePath.substring(0, dotIdx);
        return  extension;
    }

    public static File[] listFilesEndingWith(
            File dir,
            String suffix
    )
    {
        final String sfx = suffix;
        File[] tif_train = dir.listFiles(
                new FilenameFilter() {
                    public boolean accept(File dir, String name) {
                        return name.toLowerCase().endsWith(sfx);
                    }
                }
        );
        return tif_train;
    }

    public static void extractMoments2D(
            float[][]   patchin,
            double[]    center,
            double[]    theta
    )
    {

        // moments[0] = centroid(x)
        // moments[1] = centroid(y)
        // moments[2] = central_moment(x,x)
        // moments[3] = central_moment(y,y)
        // moments[4] = central_moment(x,y)

        double M00  = 0;
        double M10  = 0;
        double M01  = 0;
        double M11  = 0;
        double M20  = 0;
        double M02 	= 0;

        // 1st order... - M10, M01,...  2nd order mu11, mu20, mu02 etc...
        // image_coordinates has to be 2xN where rows 1..2 are coords
        // image_values      has to be 1xN and contains image intensities

        for (int i = 0; i < patchin.length; i++) {
            for (int j = 0; j < patchin[0].length; j++) {
                M00 += patchin[i][j];
                M10 += patchin[i][j] * i;
                M01 += patchin[i][j] * j;
                M11 += patchin[i][j] * i * j;
                M20 += patchin[i][j] * i * i;
                M02 += patchin[i][j] * j * j;
            }
        }
        // first centroids
        center[0] = M10 / M00; // centroid(x), mi10
        center[1] = M01 / M00; // centroid(y), mi01
        // second central moments
        double mu11 = M11/M00 - center[0]*center[1]; // central_moment(x,y)
        double mu20 = M20/M00 - center[0]*center[0]; // central_moment(x,x)
        double mu02 = M02/M00 - center[1]*center[1]; // central_moment(y,y)

        theta[0] = 0.5 * Math.atan((2*mu11)/(mu20-mu02));

    }

    public static float[] extractEllipse(ArrayList<int[]> locs) // row, col, majorAxis, minorAxis, angle
    {
        float[] ellipseParams = new float[5];

        double M00  = locs.size();
        double M10  = 0;
        double M01  = 0;
        double M11  = 0;
        double M20  = 0;
        double M02 	= 0;

        for (int i=0; i<locs.size(); i++) {

            M10 += locs.get(i)[0];
            M01 += locs.get(i)[1];
            M11 += locs.get(i)[1] * locs.get(i)[0];
            M20 += locs.get(i)[0] * locs.get(i)[0];
            M02 += locs.get(i)[1] * locs.get(i)[1];

        }

        ellipseParams[0] = (float) (M10 / M00);
        ellipseParams[1] = (float) (M01 / M00);

        if(M00>2) {

            float mu11 = (float) (M11 / M00) - ellipseParams[0]*ellipseParams[1];
            float mu20 = (float) (M20 / M00) - ellipseParams[0]*ellipseParams[0];
            float mu02 = (float) (M02 / M00) - ellipseParams[1]*ellipseParams[1];

            //IJ.log("EIGEN: mu20:"+mu20+",mu11:"+mu11+",mu02:"+mu02);


            ArrayList<float[]> out = eigen(mu20, mu11, mu11, mu02);

            //IJ.log("OUT");
            //for (int i=0; i<out.size(); i++) {
            //    for (int j=0; j<out.get(i).length; j++) {
            //        IJ.log(""+out.get(i)[j]+" , ");
            //    }
            //}

            if (out.size()>1) {

                // ellipseParams[2] is smaller one, ellipse[4] is it's angle
                if ( Math.abs(out.get(0)[0]) >= Math.abs(out.get(1)[0]) ) {
                    ellipseParams[2] = Math.abs(out.get(1)[0]);
                    ellipseParams[3] = Math.abs(out.get(0)[0]);
                    //IJ.log("from vec: 2_ "+out.get(1)[2]+" , 1_ "+out.get(1)[1]);
                    ellipseParams[4] = (float) (Math.atan2(out.get(1)[2], out.get(1)[1]) * (180/Math.PI));
                }
                else {
                    ellipseParams[2] = Math.abs(out.get(0)[0]);
                    ellipseParams[3] = Math.abs(out.get(1)[0]);
                    //IJ.log("from vec: 2_ "+out.get(0)[2]+" , 1_ "+out.get(0)[1]);
                    ellipseParams[4] = (float) (Math.atan2(out.get(0)[2], out.get(0)[1]) * (180/Math.PI));
                }

            }
            else {
                ellipseParams[2] = ellipseParams[3] = Math.abs(out.get(0)[0]);
                //IJ.log("from vec (both are equal): 2_ "+out.get(0)[2]+" , 1_ "+out.get(0)[1]);
                ellipseParams[4] = (float) (Math.atan2(out.get(0)[2], out.get(0)[1]) * (180/Math.PI));
            }

        }
        else {
            // there was just two of them
            ellipseParams[2] = 1;
            ellipseParams[3] = 1;
            //IJ.log("just two angle was 0 ");
            ellipseParams[4] = 0;

        }

        return ellipseParams;
    }

    public static double entropy(double[] inarray, boolean normalize)
    {
        // Shannon entropy, a measure of uncertainty
        double sum = 0;
        if(normalize) {
            for (int q = 0; q < inarray.length; q++) {
                sum += inarray[q];
            }
        }

        double entropy = 0;
        for (int q = 0; q < inarray.length; q++) {
            double p = inarray[q]/sum;
            entropy += p*Math.log(p);
        }

        return -entropy;
    }

    public static double[] kurtosis(float[] inarray, double windowRatio)
    {
        int windowSize = (int) Math.ceil(windowRatio*inarray.length);

        windowSize = (windowSize<1)? 1 : windowSize;

        double[] kurt = new double[inarray.length];

        for (int w=0; w<inarray.length; w++) {

            float m = 0;

            // mean
            for (int l=-windowSize/2; l<=windowSize/2; l++) {

                int idx = w+l;

                while(idx<0) {

                    idx += inarray.length;

                }

                while (idx>inarray.length-1) {

                    idx -= inarray.length;

                }

                m += inarray[idx];

            }

            m /= (2*(windowSize/2)+1);

            /*
            // variance
            float s2 = 0;
            for (int l=-windowSize/2; l<=windowSize/2; l++) {

                int idx = w+l;

                while(idx<0) {

                    idx += sumsPerOrt.length;

                }

                while (idx>sumsPerOrt.length-1) {

                    idx -= sumsPerOrt.length;

                }

                s2 += (sumsPerOrt[idx]-m)*(sumsPerOrt[idx]-m);

            }

            s2 /= (2*(windowSize/2)+1);
            */

            // kurtosis
            float m4 = 0;
            float m2 = 0;
            for (int l=-windowSize/2; l<=windowSize/2; l++) {

                int idx = w+l;

                while(idx<0) {

                    idx += inarray.length;

                }

                while (idx>inarray.length-1) {

                    idx -= inarray.length;

                }

                m4 += Math.pow((inarray[idx]-m), 4);
                m2 += Math.pow((inarray[idx]-m), 2);

            }

            m4 /= (2*(windowSize/2)+1);
            m2 /= (2*(windowSize/2)+1);

            kurt[w] = m4/(m2*m2) - 3;

        }

        return kurt;

    }

    public static double[] runMeanShift(
            double[] 	start,
            float[] 	inputProfile,
            int 		max_iter,
            double 	    epsilon,
            int 		h,
            double[] 	T // same length as start (this would be output)
    )
    {

        for (int i = 0; i < T.length; i++) {
            T[i] = start[i];
        }

        // performance analysis
        int countInterrupted = 0;

        for (int i = 0; i < T.length; i++) {

            int iter = 0;
            double d;

            do{

                double new_pos = runOne(T[i], h, inputProfile);

                d = Math.abs(new_pos - T[i]); // absolute diff - ok for images

//				double new_value = interp1Darray((float) new_pos, 	inputProfile);
//				double pre_value = interp1Darray((float) T[i], 		inputProfile);
//				d = Math.abs(new_value - pre_value);

                T[i] = new_pos;
                iter++;
            }
            while(iter < max_iter && d > epsilon);

            // performance analysis
            if (iter==max_iter) countInterrupted++;

        }

        System.out.println("MEAN-SHIFT: "+(countInterrupted*100f/T.length)+"% points stopped by MAX_ITER");

        return T;

    }

    public static void 	runMaxShift(
            double[] 	start,
            float[] 	inputProfile,
            int 		max_iter,
            double 	    epsilon,
            int 		h,
            double[] 	T // same length as start (this would be output)
    )
    {
        for (int i = 0; i < T.length; i++) {
            T[i] = start[i];
        }

        for (int i = 0; i < T.length; i++) {

            int iter = 0;
            double d;

            do{

                double new_pos = runOneMax(T[i], h, inputProfile);
                double new_value = interp1Darray((float) new_pos, 	inputProfile);
                double pre_value = interp1Darray((float) T[i], 		inputProfile);
                d = Math.abs(new_value - pre_value);

                T[i] = new_pos;
                iter++;
            }
            while(iter < max_iter && d > epsilon);

        }

    }

    private static double 	runOneMax(
            double  curr_pos,
            int     h,
            float[] inputProfile)
    {
        double 	    new_pos     = curr_pos;
        double 		sum 		= 0;
        double 		max	 		= Double.MIN_VALUE;

        for (float x = -h; x <= h; x+=.5f) {
            //if (true) {// ((  curr_pos + x >=0) && ( curr_pos + x <= inputProfile.length-1) ){
            //    if (true) {
            double value_read = interp1Darray((float)(curr_pos+x), inputProfile);
            if (value_read>max) {
                max = value_read;
                new_pos = curr_pos+x;
            }
            //new_pos 	+= value_read * x;
            //sum 		+= value_read 	 ;
            //}
            //}
        }
        //if(sum>0){
        //	new_pos = new_pos/sum + curr_pos;
//            new_pos[1] = new_pos[1]/sum + curr_pos[1];
        // wrap it again, this time as a position
        new_pos = (new_pos<-0.5)?(new_pos+inputProfile.length):((new_pos>=inputProfile.length-0.5)?(new_pos-inputProfile.length):new_pos);
        //}
        //else {
        //	new_pos = curr_pos;
        //}

        return new_pos;
    }

    private static double 	runOne(
            double  curr_pos,
            int     h,
            float[] inputProfile)
    {
        double 	    new_pos     = 0;
        double 		sum 		= 0;

        for (float x = -h; x <= h; x+=.5f) {
            //if (true) {// ((  curr_pos + x >=0) && ( curr_pos + x <= inputProfile.length-1) ){
            //    if (true) {
            double value_read = interp1Darray((float)(curr_pos+x), inputProfile);
            new_pos 	+= value_read * x;
            sum 		+= value_read 	 ;
            //}
            //}
        }
        if(sum>0){
            new_pos = new_pos/sum + curr_pos;
//            new_pos[1] = new_pos[1]/sum + curr_pos[1];
            // wrap it again, this time as a position
            new_pos = (new_pos<-0.5)?(new_pos+inputProfile.length):((new_pos>=inputProfile.length-0.5)?(new_pos-inputProfile.length):new_pos);
        }
        else {
            new_pos = curr_pos;
        }

        return new_pos;
    }

    public static double interp1Darray(
            float 	x,
            float[] inarray
    )
    {
        int wth = inarray.length;

        // x is a real position of a pixel in an image
        // x = [-0.5, 	W-0.5)
        // width is array width - indexes will be wrapped - first corresponds to the last one
        x = (x<=-0.5f+VERY_SMALL_POSITIVE_FLOAT		&&	x>=-0.5f-VERY_SMALL_POSITIVE_FLOAT)? -0.5f : x;
        x = (x<=wth-0.5f+VERY_SMALL_POSITIVE_FLOAT	&&	x>=wth-0.5f-VERY_SMALL_POSITIVE_FLOAT)? wth-0.5f : x;
        x = (x<-0.5f-VERY_SMALL_POSITIVE_FLOAT)? x+wth : (x>wth-0.5f+VERY_SMALL_POSITIVE_FLOAT)? x-wth : x ; // wrap

        // find four surrounding points
        int p11 = (int) Math.floor(x);
        int p12 = (int) Math.ceil(x);

        // bilinear coeffs a,b
        double a = ((p12 -p11)>0)?(x-p11)/(p12-p11) : 0.5;
        // wrap cols of surrounding pts
        //p11 = (p11<-0.5)? p11+wth : ((p11>=wth-0.5)?(p11-wth) : p11);
        p11 = (p11<0)? p11+wth : ((p11>=wth)?(p11-wth) : p11);
        p12 = (p12<0)? p12+wth : ((p12>=wth)?(p12-wth) : p12);

        double I11 = inarray[p11];
        double I12 = inarray[p12];

        return (1-a)*I11 + a*I12;
    }

    public static Vector<float[]> extractClusters(
            double[] msConvergence,
            double range,
            int M
    )
    {

        Arrays.sort(msConvergence);

        boolean[] 	seed_clustered 		= new boolean[msConvergence.length];
        int[] 		seed_cluster_idx 	= new int[msConvergence.length];
        Vector<Integer> cluster_sizes	= new Vector<Integer>();
        int nr_clusters 				= 0;

        for (int i = 0; i < msConvergence.length; i++) {
            for (int j = 0; j < msConvergence.length; j++) {

                if(!seed_clustered[j] && j!=i){

                    if(Math.abs(msConvergence[i]-msConvergence[j])<range){

                        if(seed_clustered[i]){
                            // i was clustered
                            seed_clustered[j] = true;
                            seed_cluster_idx[j] = seed_cluster_idx[i];

                            int tmp = cluster_sizes.get(seed_cluster_idx[i]);
                            tmp++;
                            cluster_sizes.setElementAt(tmp, seed_cluster_idx[i]);

                        }
                        else{
                            // i was not clustered
                            seed_clustered[i] = seed_clustered[j] = true;
                            seed_cluster_idx[i] = seed_cluster_idx[j] = nr_clusters;
                            cluster_sizes.add(2);
                            nr_clusters++;
                        }
                    }
                }
            }

            if(!seed_clustered[i]){
                seed_clustered[i] = true;
                seed_cluster_idx[i] = nr_clusters;
                cluster_sizes.add(1);
                nr_clusters++;
            }

        }

        Vector<float[]> clust = new Vector<float[]>();
        for (int k = 0; k < cluster_sizes.size(); k++) {
            if (cluster_sizes.get(k)>=M) {

                float sum = 0;
                int cnt = 0;

                for (int l=0; l<seed_cluster_idx.length; l++) {

                    if (seed_cluster_idx[l]==k) {

                        sum+=msConvergence[l];
                        cnt++;

                    }

                }

                if (cnt>0) {
                    clust.add(new float[] {sum/cnt, cnt}); // centroid, number of points
                }

            }
        }

        return clust;

    }

    public static Vector<float[]> extractClusters1(
            double[] 	msConv,
            double	    d,
            int			M,
            int         inputArrayLength // to avoid splitting in 2 clusters when values localize around the breakpoint
    )
    {

        Arrays.sort(msConv);

        boolean[] 	seed_clustered 		= new boolean[msConv.length];
        int[] 		seed_cluster_idx 	= new int[msConv.length];
        Vector<Integer> cluster_sizes	= new Vector<Integer>();
        int nr_clusters 				= 0;


        // FIRST
        seed_clustered[0] = true;
        seed_cluster_idx[0] = nr_clusters;
        cluster_sizes.add(1);
        nr_clusters++;

        // CHECK WITH LAST
        if ( Math.abs((inputArrayLength-0.5-msConv[msConv.length-1]))+Math.abs(msConv[0]+0.5) < d ) {

            seed_clustered[msConv.length-1] = true;
            seed_cluster_idx[msConv.length-1] = seed_cluster_idx[0];
            int tmp = cluster_sizes.get(seed_cluster_idx[msConv.length-1]);
            tmp++;
            cluster_sizes.setElementAt(tmp, seed_cluster_idx[msConv.length-1]);

            // correct values for calculating centroids later
//            System.out.println("DIFF "+inputArrayLength+"  btw" + msConv[msConv.length-1] + " and " + msConv[0] +" IS CALC. "+
//                    Math.abs((inputArrayLength-0.5-msConv[msConv.length-1]))
//                    +" plus "+Math.abs(msConv[0]+0.5)+" was < "+ d + " ? "+ (Math.abs((inputArrayLength-0.5-msConv[msConv.length-1]))+Math.abs(msConv[0]+0.5) < d ));
//            System.out.println("msConv "+(msConv.length-1)+"  "+msConv[msConv.length-1]+" becomes ");
            msConv[msConv.length-1] -= inputArrayLength;    // correction
//            System.out.println(" "+msConv[msConv.length-1] + " because ");

            for (int i = msConv.length-2; i>=0; i--) {
                if (Math.abs(msConv[i]-msConv[i+1])<d) {
                    seed_clustered[i] = true;
                    seed_cluster_idx[i] = seed_cluster_idx[i+1]; // TODO optimize this assignment, can be simplified!
                    tmp = cluster_sizes.get(seed_cluster_idx[i]);
                    tmp++;
                    cluster_sizes.setElementAt(tmp, seed_cluster_idx[i]);
                    msConv[i] -= inputArrayLength; // correction
                }
                else {
                    break;
                }
            }


        }

        for (int i = 1; i < msConv.length && !seed_clustered[i]; i++) {

            if(Math.abs(msConv[i]-msConv[i-1])<d){
                seed_clustered[i] = true;
                seed_cluster_idx[i] = seed_cluster_idx[i-1];
                int tmp = cluster_sizes.get(seed_cluster_idx[i]);
                tmp++;
                cluster_sizes.setElementAt(tmp, seed_cluster_idx[i]);
            }
            else {
                seed_clustered[i] = true;
                seed_cluster_idx[i] = nr_clusters;
                cluster_sizes.add(1);
                nr_clusters++;
            }

        }

        Vector<float[]> clust = new Vector<float[]>();
        for (int k = 0; k < cluster_sizes.size(); k++) {

            if (cluster_sizes.get(k)>=M) {

                float sum = 0;
                int cnt = 0;

                for (int l=0; l<seed_cluster_idx.length; l++) {

                    if (seed_cluster_idx[l]==k) {
                        //if(k==0) System.out.print(""+msConv[l]+"  , ");
                        sum+=msConv[l];
                        cnt++;

                    }

                }

                if (cnt>0) {
                    float clusterCentroid = sum/cnt;
                    clusterCentroid = (clusterCentroid<0)?                  clusterCentroid+inputArrayLength : clusterCentroid ;
                    clusterCentroid = (clusterCentroid>=inputArrayLength)?  clusterCentroid-inputArrayLength : clusterCentroid ;

//                    IJ.log("cluster index "+k+" size: "+cluster_sizes.get(k)+ " counted: " +cnt+" centroid at index "+ clusterCentroid+ " was ok with count? "+(cluster_sizes.get(k)==cnt));
                    if (!(cluster_sizes.get(k)==cnt)) {
                        System.out.println("FATAL ERROR!");
                    }
//                    System.out.println("add "+clusterCentroid+" at pos "+clust.size());
                    clust.add(new float[] {clusterCentroid, cnt}); // centroid, number of points
                }

            }
        }

        return clust;

    }

    /**
     * background median estimation ( Wirth's algorithm )
     * Title: Algorithms + data structures = programs
     * Publisher: Englewood Cliffs: Prentice-Hall, 1976
     * Physical description: 366 p.
     * Series: Prentice-Hall Series in Automatic Computation
     */
    public static double median_Wirth(float[] a)
    {
        int n = a.length;
        int i, j, l, m, k;
        double x;
        if (n % 2 == 0) k = (n/2)-1;
        else k = (n/2);
        l=0 ; m=n-1 ;
        while (l < m)
        {
            x=a[k] ;
            i = l ;
            j = m ;
            do
            {
                while (a[i] < x) i++ ;
                while (x < a[j]) j-- ;
                if (i <= j) {
                    float temp = a[i];
                    a[i] = a[j];
                    a[j] = temp;
                    i++ ; j-- ;
                }
            } while (i <= j) ;
            if (j < k) l = i ;
            if (k < i) m = j ;
        }
        return a[k] ;
    }

    public static float[] hdome(float[] I, float h)
    {

        // sequential implementation (Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms. Vincent 1993.)
        float[] J       = new float[I.length];
        float[] Jprev   = new float[I.length];

        // create marker image
        for (int w=0; w<I.length; w++)
            J[w]        = I[w] - h;

        float diff;

        do {

            for (int i1=0; i1<J.length; i1++)
                Jprev[i1] = J[i1];

            // calculate next J
            // raster order
            for (int i1=0; i1<J.length; i1++)
                J[i1] = (i1>0)? Math.min(Math.max(J[i1-1], J[i1]), I[i1]) : Math.min(J[i1], I[i1]) ;
            // anti-raster
            for (int i1=J.length-1; i1>=0; i1--)
                J[i1] = (i1<J.length-1)? Math.min(Math.max(J[i1+1], J[i1]), I[i1]) : Math.min(J[i1], I[i1]) ;

            // calculate diff
            diff = 0;
            for (int i1=0; i1<J.length; i1++) {
                diff += Math.abs(J[i1]-Jprev[i1]);
            }

        }
        while (diff>1);

        // subtract the reconstruction
        for (int i1=0; i1<I.length; i1++) {
            J[i1] = I[i1] - J[i1];
        }

        return J;
    }

    public static float[] hdome_Circular(float[] I, float h, float hmin)
    {

        // sequential implementation (Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms. Vincent 1993.)
        float[] J       = new float[I.length];
        float[] Jprev   = new float[I.length];

        // create marker image
        for (int w=0; w<I.length; w++)
            J[w]        = I[w] - h;

        float diff;

        do {

            for (int i1=0; i1<J.length; i1++)
                Jprev[i1] = J[i1];

            // calculate next J
            // raster order
            for (int i1=0; i1<J.length; i1++)
                J[i1] = (i1>0)? Math.min(Math.max(J[i1-1], J[i1]), I[i1]) : Math.min(Math.max(J[J.length-1], J[i1]), I[i1]) ;
            // anti-raster
            for (int i1=J.length-1; i1>=0; i1--)
                J[i1] = (i1<J.length-1)? Math.min(Math.max(J[i1+1], J[i1]), I[i1]) : Math.min(Math.max(J[0], J[i1]), I[i1]) ;

            // calculate diff
            diff = 0;
            for (int i1=0; i1<J.length; i1++) {
                diff += Math.abs(J[i1]-Jprev[i1]);
            }

        }
        while (diff>1);

        // subtract the reconstruction
        for (int i1=0; i1<I.length; i1++) {
            J[i1] = (I[i1] - J[i1]>hmin)? 1 : 0 ;
        }

//        for (int i1=0; i1<I.length; i1++) {
//            J[i1] = (J[i1]>hmin)? 1 : 0;
//        }

        return J;
    }

    public static float wrap_180(float angDeg) // -180 + 180 wrapping
    {

        float out = angDeg;

        while (out<-180) {
            out += 360;
        }

        while (out>=180) {
            out -= 360;
        }

        return out;

    }

    public static float wrap_360(float angDeg)
    {
        float out = angDeg;

        while (out<0) {
            out += 360;
        }

        while (out>=360) {
            out -= 360;
        }

        return out;
    }

    public static double min3(double a, double b, double c)
    {
        return Math.min(a, Math.min(b, c));
    }

    public static float min3(float[] a)
    {
        return Math.min(a[0], Math.min(a[1], a[2]));
    }

    public static float max3(float[] a)
    {
        return Math.max(a[0], Math.max(a[1], a[2]));
    }

    public static PolygonRoi drawEllipse(float x, float y, float a, float b, float angle, Color clr, float lineWidth, int nrPoints)
    {

        float[] xEll = new float[nrPoints];
        float[] yEll = new float[nrPoints];
        double step = (2*Math.PI)/nrPoints;
        double beta = -angle*(Math.PI/180);

        for (int i=0; i<nrPoints; i++) {
            double alpha = i*step;
            xEll[i] = (float) (x + a * Math.cos(alpha) * Math.cos(beta) - b * Math.sin(alpha) * Math.sin(beta));
            yEll[i] = (float) (y + a * Math.cos(alpha) * Math.sin(beta) + b * Math.sin(alpha) * Math.cos(beta));
        }

        PolygonRoi p = new PolygonRoi(xEll, yEll, nrPoints, Roi.POLYGON);
        p.setStrokeWidth(lineWidth);
        p.setStrokeColor(clr);

        return p;

    }

    public static ArrayList<float[]> eigen(float a11, float a12, float a21, float a22)
    {

        float a = 1;
        float b = -a11-a22;
        float c = a11*a22 - a12*a21;

        double discriminant = b*b-4*a*c;

        ArrayList<float[]> out = new ArrayList<float[]>();

        if (discriminant>VERY_SMALL_POSITIVE) {

            float norml, v1_lmb1, v2_lmb1, v1_lmb2, v2_lmb2;

            // 2 distinct real roots
            float lambda1 = (float) ((-b + Math.sqrt(discriminant)) / (2*a));
            float lambda2 = (float) ((-b - Math.sqrt(discriminant)) / (2*a));

            if (a12<VERY_SMALL_POSITIVE && a12 >-VERY_SMALL_POSITIVE) {

                //a12~0
                v1_lmb1 = 0;
                v2_lmb1 = 1;

                v1_lmb2 = -1;
                v2_lmb2 = 0;

            }
            else {

                v1_lmb1 = 1;
                v2_lmb1 = (float) ((lambda1-a11)/a12);
                // normalize them
                norml = (float) Math.sqrt(1+v2_lmb1*v2_lmb1);
                v1_lmb1 = v1_lmb1 / norml;
                v2_lmb1 = v2_lmb1 / norml;

                // lambda2 vectors
                v1_lmb2 = 1;
                v2_lmb2 = (float) ((lambda2-a11)/a12);
                // normalize them
                norml = (float) Math.sqrt(1+v2_lmb2*v2_lmb2);
                v1_lmb2 = v1_lmb2 / norml;
                v2_lmb2 = v2_lmb2 / norml;

            }

            out.add(new float[]{lambda1, v1_lmb1, v2_lmb1});
            out.add(new float[]{lambda2, v1_lmb2, v2_lmb2});

        }
        else if (discriminant<-VERY_SMALL_POSITIVE) {
            // complex roots
            out.add(new float[]{0, 0, 0});

        }
        else {
            // one real root

            float norml, v1_lmb1, v2_lmb1, v1_lmb2, v2_lmb2;

            float lambda1 = (float) ((-b) / (2*a));

            if (a12<VERY_SMALL_POSITIVE && a12 >-VERY_SMALL_POSITIVE) {

                //a12~0
                v1_lmb1 = 0;
                v2_lmb1 = 1;

                v1_lmb2 = -1;
                v2_lmb2 = 0;

            }
            else {

                // lambda1 vectors
                v1_lmb1 = 1;
                v2_lmb1 = (float) ((lambda1-a11)/a12);
                // normalize them
                norml = (float) Math.sqrt(1+v2_lmb1*v2_lmb1);
                v1_lmb1 = v1_lmb1 / norml;
                v2_lmb1 = v2_lmb1 / norml;



            }



            out.add(new float[]{lambda1, v1_lmb1, v2_lmb1});

        }

        return out;

    }

    public static ArrayList<int[]> comb(int n, int k)
    {

        ArrayList<int[]> out = new ArrayList<int[]>();
        int[] cmbs = new int[k];

        // initial combination
        for (int i=0; i<k; i++) {
            cmbs[i] = i;
        }

        // add initial combination
        int[] cmbToAdd;

        cmbToAdd = new int[k];
        for (int i=0; i<k; i++) cmbToAdd[i] = cmbs[i];
        out.add(cmbToAdd);

        while (nextComb(cmbs, n, k)) {
            cmbToAdd = new int[k];
            for (int i=0; i<k; i++) cmbToAdd[i] = cmbs[i];
            out.add(cmbToAdd);
        }

        return out;

    }

    public static float[] getMinMax(float[] a)
    {

        float[] minMax = new float[2];
        minMax[0] = Float.MAX_VALUE;
        minMax[1] = Float.MIN_VALUE;
        for (int i = 0; i < a.length; i++) {
            if (a[i]<minMax[0]) minMax[0] = a[i];
            if (a[i]>minMax[1]) minMax[1] = a[i];
        }
        return minMax;

    }

    public static float findNextStationaryValue(float pos, int[] roundedPeaks, float[] profile) // , float sidePos1, float sidePos2,
    {
        // 2 threads, each initializes when the next one is lower
        int l = Tools.wrap( (int) Math.floor(pos), profile.length);
        int r = Tools.wrap( (int) Math.ceil(pos), profile.length);

        // align
        float diff;
        diff = profile[l]-profile[Tools.wrap(l+1, profile.length)];
        while (diff>=0) {
            l--;
            l = Tools.wrap(l, profile.length);
            diff = profile[l]-profile[Tools.wrap(l+1, profile.length)];
        }
        diff = profile[r]-profile[Tools.wrap(r-1, profile.length)];
        while (diff>=0) {
            r++;
            r = Tools.wrap(r, profile.length);
            diff = profile[r] - profile[Tools.wrap(r-1, profile.length)];
        }

        // step ahead
        l--;
        l = Tools.wrap(l, profile.length);
        r++;
        r = Tools.wrap(r, profile.length);

        diff = profile[l]-profile[Tools.wrap(l+1, profile.length)];
        while (diff<0) {

            // check
            for (int g=0; g<roundedPeaks.length; g++) if (l==roundedPeaks[g]) return profile[Math.round(pos)];  // in case it reaches other 2 peaks

            l--;
            l = Tools.wrap(l, profile.length);
            diff = profile[l]-profile[Tools.wrap(l+1, profile.length)];
        }

        diff = profile[r]-profile[Tools.wrap(r-1, profile.length)];
        while (diff<0) {

            // check
            for (int g=0; g<roundedPeaks.length; g++) if (r==roundedPeaks[g]) return profile[Math.round(pos)];  // in case it reaches other 2 peaks

            r++;
            r = Tools.wrap(r, profile.length);
            diff = profile[r]-profile[Tools.wrap(r-1, profile.length)];
        }

        return (profile[l]>profile[r])?
                Math.min(profile[Math.round(pos)], profile[l]) :
                Math.min(profile[Math.round(pos)], profile[r]) ;
        //((detection[l]<detection[Math.round(pos)])?detection[l]:0)     :
        //((detection[r]<detection[Math.round(pos)])?detection[r]:0)     ; // take the one that is closer to the detection[pos]

    }

    public static int wrap(int idx, int len)
    {
        return idx = (idx<0)? idx+len : (idx>=len)? idx-len : idx ;
    }

    public static float[] circularLocalAvg(float[] in, int w)
    {
        float[] out = new float[in.length];
        float locAvg;
        int wrappedIdx;
        for (int i=0; i<in.length; i++) {
            locAvg = 0;
            for (int j=i-w; j<=i+w; j++) {
                wrappedIdx = (j<0)? j+in.length : (j>=in.length)? j-in.length : j ;
                locAvg += in[wrappedIdx];
            }
            out[i] = locAvg/(2*w+1);
        }
        return out;
    }

    public static int[][] hungarianMappingAnglesDeg(ArrayList<Float> angsDegA, ArrayList<Float> angsDegB)  // will output the best matching pairs of indexes
    {

        // match A and B dimentional arrays using Euclidean distance, output the mapping as array of correspondence indexes
        boolean[][]     chkd = new boolean[angsDegA.size()][angsDegB.size()]; // marker
        double[][]      dst2 = new double[angsDegA.size()][angsDegB.size()];

        for (int i = 0; i < angsDegA.size(); i++) {
            for (int j = 0; j < angsDegB.size(); j++) {
                dst2[i][j] = Math.abs(Tools.wrap_180(angsDegA.get(i) - angsDegB.get(j)));
            }
        }

        int totalPairs = Math.min(angsDegA.size(), angsDegB.size());

        int[][] pairs = new int[totalPairs][2];

        for (int check=0; check<totalPairs; check++) {

            double dst2Min = Double.MAX_VALUE;
            int imin = -1;
            int jmin = -1;

            for (int i=0; i<angsDegA.size(); i++) {
                for (int j=0; j<angsDegB.size(); j++) {
                    if (!chkd[i][j] && dst2[i][j]<dst2Min) {
                        dst2Min = dst2[i][j];
                        imin = i;
                        jmin = j;
                    }
                }

            }

            // row imin in chkd to true
            for (int w=0; w<angsDegB.size(); w++) chkd[imin][w] = true;  // loop columns
            // col jmin in chkd to true
            for (int w=0; w<angsDegA.size(); w++) chkd[w][jmin] = true;  // loop rows

            pairs[check][0] = imin;
            pairs[check][1] = jmin;

//            IJ.log("matched angle "+angsDegA.get(imin)+" and "+angsDegB.get(jmin)+" diff="+dst2Min);

        }

        return pairs;

    }

    public static void swap(float[] vals, int[] mapping)
    {

        float[] temp = new float[vals.length];

        for (int t = 0; t < vals.length; t++) {
            temp[mapping[t]] = vals[t];
        }

        for (int t =0; t < vals.length; t++) {
            vals[t] = temp[t];
        }

    }

    public static int[] hungarian33(float[] refAng, float[] currAng)
    {
        //int[][] mapping = new int[3][3];

        // match using Hungarian algorithm
        boolean[][] chkd = new boolean[3][3];
        double[][] 	dst2 = new double[3][3];

        for (int i=0; i<3; i++) {
            for (int j=0; j<3; j++) {
                dst2[i][j] = Math.abs(Tools.wrap_180(currAng[i] - refAng[j]));
            }
        }

        int[] mapping = new int[3];

        for (int check=0; check<3; check++) {

            double dst2Min = Double.MAX_VALUE;
            int imin = -1;
            int jmin = -1;

            for (int i=0; i<3; i++) {

                for (int j=0; j<3; j++) {
                    if (!chkd[i][j] && dst2[i][j]<dst2Min) {
                        dst2Min = dst2[i][j];
                        imin = i;
                        jmin = j;
                    }
                }

            }

            // row imin in chkd to true
            for (int w=0; w<3; w++) chkd[imin][w] = true;
            // col jmin in chkd to true
            for (int w=0; w<3; w++) chkd[w][jmin] = true;

            mapping[imin]=jmin;

        }

        return mapping;

    }

    public static void swap3(float[] vals, int[] mapping)
    {

        float[] temp   = new float[3]; // ugly solution to swap 3
        //float[] indexesTemp     = new float[3]; // ugly solution to swap 3
        //float[] peaksTemp       = new float[3];

        for (int t=0; t<3; t++) {
            temp[mapping[t]] = vals[t];
//            indexesTemp[mapping[t]] = indexes[t];
//            peaksTemp[mapping[t]]   = peaks[t];
        }

//                        IJ.log("...");
//                        IJ.log(Arrays.toString(anglesDegTemp));
//                        IJ.log("...");

        for (int t=0; t<3; t++) {
            vals[t]      = temp[t];
//            anglesDeg[t]    = anglesDegTemp[t];
//            peaks[t]        = peaksTemp[t];
        }

    }

    private static boolean nextComb(int[] cmbs, int n, int k)
    {

        int i = k-1;
        ++cmbs[i];
        while (i>0 && cmbs[i] >= n-k+1+i) {
            --i;
            ++cmbs[i];
        }

        if(cmbs[0]>n-k) {
            return false;
        }

        for (i=i+1; i<k; i++) {
            cmbs[i] = cmbs[i-1] + 1;
        }

        return true;
    }

    public static float triangleArea(float xA, float yA, float xB, float yB, float xC, float yC)
    {
        return (float) (0.5*Math.abs((xA-xC)*(yB-yA)-(xA-xB)*(yC-yA)));
    }

    public static float angularDeviation(float Cx, float Cy, float Ax, float Ay, float Bx, float By)
    {
        float dotProd = 0;
        float v12 = 0;
        float v32 = 0;
        dotProd = (Ax-Cx) * (Bx-Ax) + (Ay-Cy) * (By-Ay);
        v12 = (float) Math.sqrt((Ax-Cx)*(Ax-Cx) + (Ay-Cy)*(Ay-Cy));
        v32 = (float) Math.sqrt((Bx-Ax)*(Bx-Ax) + (By-Ay)*(By-Ay));
        return dotProd/(v12*v32);
    }

    public static float angularDeviation(float[] x1, float[] x2, float[] x3)
    {
        float dotProd = 0;
        float v12 = 0;
        float v32 = 0;
        for (int i = 0; i < x1.length; i++) {
            dotProd += (x2[i]-x1[i]) * (x3[i]-x2[i]);
            v12 += (x2[i]-x1[i]) * (x2[i]-x1[i]);
            v32 += (x3[i]-x2[i]) * (x3[i]-x2[i]);
        }
        return (float) (dotProd/(Math.sqrt(v12)*Math.sqrt(v32)));
    }

    public static String createDir(String dir_name)
    {

        File f = new File(dir_name);

        if(!f.exists()){  //

            try{
                // Create one directory
                boolean success = f.mkdirs();
                if (!success) {
                    System.out.println("Error: Directory: " + dir_name + "    .... was NOT created");
                }

            }
            catch (Exception e){//Catch exception if any
                System.err.println("Error: " + e.getMessage());
            }

        }
        else {



            // empty it before deleting
            String[] entries = f.list();

            for(String s: entries){
                File currentFile = new File(f.getPath(),s);
                currentFile.delete();
            }

            f.delete();

            try{
                // Create one directory
                boolean success = f.mkdirs();
                if (!success) {
                    System.out.println("Error: Directory: " + dir_name + "    .... was NOT created");
                }

            }
            catch (Exception e){//Catch exception if any
                System.err.println("Error: " + e.getMessage());
            }

        }

        return f.getAbsolutePath();

    }

    public static void createDirs(String dirs_names)
    {
        try{
            // Create multiple directories
            boolean success = (new File(dirs_names)).mkdirs();
            if (success) {
                System.out.println("Directories: " + dirs_names + " created");
            }
        }
        catch (Exception e){//Catch exception if any
            System.err.println("Error: " + e.getMessage());
        }
    }

    public static int[] descending(ArrayList<Integer> a)
    {

        // prepare array with indexes first
        int[] idx = new int[a.size()];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.size()-1; i++) {
            for (int j = i+1; j < a.size(); j++) {
                if (a.get(j)>a.get(i)) { // desc.
                    int temp 	= a.get(i);
                    a.set(i, a.get(j));
                    a.set(j, temp);

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

    public static int[] descending(float[] a)
    {

        // prepare array with indexes first
        int[] idx = new int[a.length];
        for (int i=0; i<idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.length-1; i++) {
            for (int j = i+1; j < a.length; j++) {
                if (a[j]>a[i]) { // desc.
                    float temp 	= a[i];
                    a[i]		= a[j];
                    a[j] 		= temp;

                    int temp_idx 	= idx[i];
                    idx[i] 			= idx[j];
                    idx[j]			= temp_idx;
                }
            }
        }

        return idx;

    }

    public static int[] descendingFloat(ArrayList<Float> a)
    {

        int[] idx = new int[a.size()];

        for (int i = 0; i < idx.length; i++) idx[i] = i;

        for (int i = 0; i < a.size()-1; i++) {
            for (int j = i+1; j < a.size(); j++) {
                if (a.get(j)>a.get(i)) {
                    float temp = a.get(i);
                    a.set(i, a.get(j));
                    a.set(j, temp);

                    int tmp_idx = idx[i];
                    idx[i]      = idx[j];
                    idx[j]      = tmp_idx;
                }
            }
        }

        return idx;

    }

    public static String getFileExtension(String file_path)
    {
        String extension = "";

        int i = file_path.lastIndexOf('.');
        if (i >= 0) {
            extension = file_path.substring(i+1);
        }

        return extension;
    }

    public static String getFileName(String file_path)
    {
        String name = "";

        int i = file_path.lastIndexOf('.');
        int j = file_path.lastIndexOf(File.separator);

        if (i > 0) {
            if (j>=0) {
                name = file_path.substring(j+1, i);
            }
            else {
                name = file_path.substring(0, i);
            }
        }

        return name;
    }
}
