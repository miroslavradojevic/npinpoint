package com.braincadet.npinpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Prefs;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

/** */
public class ProfileSpread extends Thread {

    private int begN, endN;

    // input
    public static short[][]	    prof2;                      // profiles
    public static int[][] 	    i2xy;                       // index to xy location
    public static int[][]     	xy2i;                       // xy location to index
    public static int W;
    public static int H;

    // output - profile spread in some forms
    public static float[] 		spread;                     // measure of spread of profile here (variance or range)
    public static boolean[] 	profile_diverse;			// profile has enough spread

    public ProfileSpread(int n0, int n1) {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(int[][] _i2xy, int[][] _xy2i, short[][] _prof2, int _W, int _H){

        i2xy      		= _i2xy;
        prof2           = _prof2;
        xy2i 			= _xy2i;
        W = _W;
        H = _H;

        // allocate output
        spread = new float[i2xy.length];
        profile_diverse = new boolean[i2xy.length];

    }

    public void run() {

        int len = prof2[0].length;
        float[] f = new float[len]; // to store the profiles

        //main
        for (int locationIdx = begN; locationIdx < endN; locationIdx++) { // all foreground locations


            float profile_max = Float.NEGATIVE_INFINITY;
            float profile_min = Float.POSITIVE_INFINITY;

            for (int i=0; i<len; i++) {

                f[i] = ((prof2[locationIdx][i] & 0xffff) / 65535f) * 255f; // retrieve the profile

                if (f[i]>profile_max) profile_max = f[i];
                if (f[i]<profile_min) profile_min = f[i];

            }

            spread[locationIdx] = profile_max - profile_min; //(float) Math.pow(Stat.std(f), 2);

        }

    }

    public static void threshold(){



        // will threshold the spread FloatProcessor using available automatic threshold algorithm
        ImagePlus  spread_ip = new ImagePlus("", getSpread());
        IJ.setAutoThreshold(spread_ip, "Percentile dark");
        Prefs.blackBackground = true;
        IJ.run(spread_ip, "Convert to Mask", "");
        // output is in spread_ip


        // remove those spreads that are having low variation
        // it is a matter of thresholding the histogram so that it separates two blocks well


        // store the thresholded values back to the array
        for (int i = 0; i < profile_diverse.length; i++) {
            int x = i2xy[i][0];
            int y = i2xy[i][1];

            if (spread_ip.getProcessor().get(x, y)>0) {
                profile_diverse[i] = true;
            }
            else {
                profile_diverse[i] = false;
            }


        }

    }

    public static ByteProcessor getMask() {

        ByteProcessor bp = new ByteProcessor(W, H);
        for (int i = 0; i < profile_diverse.length; i++) {

            int x = i2xy[i][0];
            int y = i2xy[i][1];

            bp.set(x, y, profile_diverse[i]?255:0);

        }

        return bp;

    }

    public static ImageStack getSpreadDistribution(int nr_bins){

        ImageStack is_out = new ImageStack(528, 255);

        float[] bins = new float[nr_bins];
        float[] distribution = new float[nr_bins];

        Hist.getDistribution(spread, nr_bins, bins, distribution);
        Plot p = new Plot("", "", "", bins, distribution, Plot.LINE);

        is_out.addSlice("", p.getProcessor());
        return is_out;

    }

    public static ImageProcessor getSpread() {
        FloatProcessor fp = new FloatProcessor(W, H);

        for (int i = 0; i < spread.length; i++) {

            int atX = i2xy[i][0];
            int atY = i2xy[i][1];

            fp.setf(atX, atY, spread[i]);

        }

        return  fp;
    }

    public static int getNrCritpointCandidates(){
        int count = 0;
        for (int i = 0; i < profile_diverse.length; i++) {
            if (profile_diverse[i])
                count++;
        }
        return count;
    }

}
