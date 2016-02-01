package com.braincadet.npinpoint;

import java.util.Arrays;

/** collection of statistical functions used in detection */
public class Stat {

    public static float median(float[] a)
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

    public static float quantile(float[] a, int ratioNum, int ratioDen) // ratioNum/ratioDen first ones
    {
        int n = a.length;
        int i, j, l, m, k;
        double x;

        if ((ratioNum*n) % ratioDen == 0) k = ((ratioNum*n)/ratioDen)-1;
        else k = (ratioNum*n)/ratioDen;

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


    public static float average(float[] a)
    {
        float meanVal = 0;
        for (int i=0; i<a.length; i++)
            meanVal += a[i];
        return meanVal/a.length;

    }

    public static float std(float[] in, float avg)
    {
        float std = 0;
        for (int i=0; i<in.length; i++) {
            std += (in[i]-avg)*(in[i]-avg);
        }
        std /= in.length;
        std = (float) Math.sqrt(std);
        return std;
    }

    public static float std(float[] in)
    {

        float avg = 0;
        for (int i=0; i<in.length; i++) {
            avg += in[i];
        }
        avg /= in.length;

        float std = 0;
        for (int i=0; i<in.length; i++) {
            std += (in[i]-avg)*(in[i]-avg);
        }
        std /= in.length;
        std = (float) Math.sqrt(std);
        return std;
    }

    public static float var(float[] in, float avg)
    {
        float var = 0;
        for (int i=0; i<in.length; i++) {
            var += (in[i]-avg)*(in[i]-avg);
        }
        var /= in.length;
        return var;
    }

    public static float var(float[] in)
    {

        float avg = 0;
        for (int i=0; i<in.length; i++) {
            avg += in[i];
        }
        avg /= in.length;

        float var = 0;
        for (int i=0; i<in.length; i++) {
            var += (in[i]-avg)*(in[i]-avg);
        }
        var /= in.length;
        return var;
    }

    public static final float get_max(float[] in) {
        float curr_max = in[0];
        for (int i=1; i<in.length; i++) {
            if (in[i]>curr_max) curr_max = in[i];
        }
        return curr_max;
    }

    public static final float get_min(float in[]) {
        float curr_min = in[0];
        for (int i=1; i<in.length; i++) {
            if (in[i]<curr_min) curr_min = in[i];
        }
        return curr_min;
    }

    public static final void min_max_normalize(float in[]) {
        float curr_min = in[0];
        float curr_max = in[0];
        for (int i=1; i<in.length; i++) {
            if (in[i]>curr_max) curr_max = in[i];
            if (in[i]<curr_min) curr_min = in[i];
        }
        if (Math.abs(curr_max-curr_min)<=Float.MIN_VALUE) {
            Arrays.fill(in, 0); // special case
        }
        for (int i=0; i<in.length; i++) {
            in[i] = (in[i]-curr_min)/(curr_max-curr_min);
            in[i] = (in[i]>1)? 1 : (in[i]<0)? 0 : in[i];
        }


    }

    public static final float[] min_max_normalize1(float[] in)
    {

        float[] out = new float[in.length];

        float curr_min = in[0];
        float curr_max = in[0];
        out[0] = in[0];

        for (int i=1; i<in.length; i++) {
            out[i] = in[i];
            if (in[i]>curr_max) curr_max = in[i];
            if (in[i]<curr_min) curr_min = in[i];
        }
        if (Math.abs(curr_max-curr_min)>Float.MIN_VALUE) {

            for (int i=0; i<out.length; i++) {
                out[i] = (out[i]-curr_min)/(curr_max-curr_min);
                out[i] = (out[i]>1)? 1 : (out[i]<0)? 0 : out[i];
            }
        }

        return out;

    }

    public static final void probability_distribution(float[] in) {

        float sum = 0;
        for (int i = 0; i < in.length; i++) {
            sum += in[i];
        }

        if (sum==0) {
            for (int i = 0; i < in.length; i++) {
                in[i] = 1f / in.length; // they are all equal weight input has all zeros
            }
        }
        else {
            for (int i = 0; i < in.length; i++) {
                in[i] /= sum;
            }
        }

    }

    public static final void normalize(float in[]) { // will normalize with respect to the max.

        float curr_max = in[0];

        for (int i = 1; i < in.length; i++) {
            if(in[i]>curr_max) curr_max = in[i];
        }

        if (curr_max<=Float.MIN_VALUE) {
            Arrays.fill(in, 0);
        }

        for (int i=0; i<in.length; i++) {
            in[i] = in[i]/curr_max;
        }

    }

    public static final float[] get_min_max(float[] in) {
        float[] out_min_max = new float[2];//min->0, max->1
        out_min_max[0] = in[0];
        out_min_max[1] = in[0];
        for (int i=1; i<in.length; i++) {
            if (in[i]<out_min_max[0]) {
                out_min_max[0] = in[i];
            }
            if (in[i]>out_min_max[1]) {
                out_min_max[1] = in[i];
            }
        }
        return out_min_max;
    }

    public static final float cemd(float f1, float f2, float f3, float f4, float g1, float g2, float g3, float g4){


        // cummulative
        float[][] F = new float[][]{
                {f1,            f1+f2,          f1+f2+f3,       f1+f2+f3+f4},
                {f1+f2+f3+f4,   f2,             f2+f3,          f2+f3+f4},
                {f1+f3+f4,      f1+f2+f3+f4,    f3,             f3+f4},
                {f1+f4,         f1+f2+f4,       f1+f2+f3+f4,    f4}
        };

        float[][] G = new float[][]{
                {g1,            g1+g2,          g1+g2+g3,       g1+g2+g3+g4},
                {g1+g2+g3+g4,   g2,             g2+g3,          g2+g3+g4},
                {g1+g3+g4,      g1+g2+g3+g4,    g3,             g3+g4},
                {g1+g4,         g1+g2+g4,       g1+g2+g3+g4,    g4}
        };

        float out = Float.POSITIVE_INFINITY;

        for (int i=0; i<F.length; i++) {

            float avg = 0;
            for (int j=0; j<F[i].length; j++) {
                avg += Math.abs(F[i][j]-G[i][j]);
            }
            avg /= F[i].length;

            if (avg<out) out = avg;

        }

        return out;

    }

    private static void extract_cemd(float f1, float f2, float f3, float f4, float[] cemd_end_bdy_bif_crs) {

        // store 4 distances
        cemd_end_bdy_bif_crs[0] = Stat.cemd(f1, f2, f3, f4, 1f, 0, 0, 0);
        cemd_end_bdy_bif_crs[1] = Math.min(Stat.cemd(f1, f2, f3, f4, 0.5f, 0.5f, 0, 0), Stat.cemd(f1, f2, f3, f4, 0.5f, 0, 0.5f, 0));
        cemd_end_bdy_bif_crs[2] = Stat.cemd(f1, f2, f3, f4, 1f/3, 1f/3, 1f/3, 0);
        cemd_end_bdy_bif_crs[3] = Stat.cemd(f1, f2, f3, f4, 1f/4, 1f/4, 1f/4, 1f/4);
    }

}
