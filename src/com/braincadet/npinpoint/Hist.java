package com.braincadet.npinpoint;

import ij.gui.Plot;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/** histogram tools */
public class Hist {

    // java -cp "$HOME/critpoint/*:$HOME/jarlib/*" Hist
    public static void main(String[] args) {

        System.out.println("test histogram...");

        int N = 10000;
        float range = 100;
        float mean = 50;
        float std = 5;

        float[] values = new float[N];

        Random rd = new Random();
        for (int i = 0; i < values.length; i++) {
            //values[i] = rd.nextFloat() * (range);
            values[i] = (float) (rd.nextGaussian() * std + mean);

        }

        // plot values
        float[] xaxis_values = new float[values.length];
        for (int i = 0; i < xaxis_values.length; i++) {
            xaxis_values[i] = i;
        }
        Plot p = new Plot("values", "", "", xaxis_values, values, Plot.LINE);
        p.show();


        // plot histogram counts
        int nr_bins = 15;

        float[] 	binx 	= new float[nr_bins];
        int[] 		count 	= new int[nr_bins];
        float[] 	distr 	= new float[nr_bins];

        Arrays.fill(binx, 0);

        getCounts(values, nr_bins, binx, count);
        System.out.println(Arrays.toString(binx));
        System.out.println(Arrays.toString(count));


        getDistribution(values, nr_bins, binx, distr);
        System.out.println(Arrays.toString(distr));

        Plot h = new Plot("histogram", "", "");
        h.setLimits(binx[0], binx[binx.length-1], 0, find_max(distr));
        h.draw();
        h.setLineWidth(3);
        h.addPoints(binx, distr, Plot.LINE);

        h.draw();
        h.show();


    }

    public static void getCounts(float[] vals, int nr_bins, float[] binx, int[] biny) {

        float max_val = find_max(vals);
        float min_val = find_min(vals);

        float step = (max_val - min_val) / nr_bins;

        for (int i = 0; i < nr_bins; i++) binx[i] = min_val + step / 2f + i * step;

        Arrays.fill(biny, 0);
        for (int i = 0; i < vals.length; i++) {
            int where = (int) Math.floor( (vals[i]-min_val) / step);
            if (where>=nr_bins) where = nr_bins - 1;
            biny[where]++;
        }

    }

    public static void getDistribution(float[] vals, int nr_bins, float[] binx, float[] distr) {

        float max_val = find_max(vals);
        float min_val = find_min(vals);

        float step = (max_val - min_val) / nr_bins;

        for (int i = 0; i < nr_bins; i++) binx[i] = min_val + step / 2f + i * step;

        Arrays.fill(distr, 0);
        int cnt = 0;
        for (int i = 0; i < vals.length; i++) {
            int where = (int) Math.floor( (vals[i]-min_val) / step);
            if (where>=nr_bins) where = nr_bins - 1;
            distr[where]++;
            cnt++;
        }

        // normalize
        for (int i = 0; i < distr.length; i++) {
            distr[i] = distr[i] / cnt;
        }

    }

    public static void getDistribution(ArrayList<Float> vals, int nr_bins, float[] binx, float[] distr) {

        float max_val = find_max(vals);
        float min_val = find_min(vals);

        float step = (max_val - min_val) / nr_bins;

        for (int i = 0; i < nr_bins; i++) binx[i] = min_val + step / 2f + i * step;

        Arrays.fill(distr, 0);
        int cnt = 0;
        for (int i = 0; i < vals.size(); i++) {
            int where = (int) Math.floor( (vals.get(i)-min_val) / step);
            if (where>=nr_bins) where = nr_bins - 1;
            distr[where]++;
            cnt++;
        }

        // normalize
        for (int i = 0; i < distr.length; i++) {
            distr[i] = distr[i] / cnt;
        }

    }

    public static float find_max(float[] values) {

        float max_value = values[0];

        for (int i = 1; i < values.length; i++) {

            if (values[i]>max_value) {
                max_value = values[i];
            }

        }

        return max_value;

    }

    public static float find_max(ArrayList<Float> values) {

        float max_value = values.get(0);

        for (int i = 1; i < values.size(); i++) {

            if (values.get(i) > max_value) {
                max_value = values.get(i);
            }

        }

        return max_value;

    }

    public static float find_min(float[] values) {

        float min_value = values[0];

        for (int i = 1; i < values.length; i++) {

            if (values[i] < min_value) {
                min_value = values[i];
            }

        }

        return min_value;

    }

    public static float find_min(ArrayList<Float> values) {

        float min_value = values.get(0);

        for (int i = 1; i < values.size(); i++) {

            if (values.get(i) < min_value) {
                min_value = values.get(i);
            }

        }

        return min_value;

    }

}
