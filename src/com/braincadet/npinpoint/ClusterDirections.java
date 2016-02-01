package com.braincadet.npinpoint;

import java.util.ArrayList;
/** cluster directions using mean-shift */
public class ClusterDirections {

    public static ArrayList<float[]> run(ArrayList<float[]> vxy, float alfaDeg) {

        int max_iter = 20;
        float epsilon = 0.001f;
        int min_cluster_cnt = 1; // to avoid taking into account outliers
        float alfaRad = Deg2Rad(alfaDeg);

        float[][] vxy_conv = meanShift(vxy, max_iter, epsilon, alfaRad); 				    // mean-shift, vxy are normalized
        int[] out_lab = clustering(vxy_conv, alfaRad); 									    // cluster convergences
        ArrayList<float[]> clustered_vxy_cnt = extracting(out_lab, vxy, min_cluster_cnt);   // extract clusters (sorted by size)

        return clustered_vxy_cnt;

    }

    private static float[][]	meanShift(ArrayList<float[]> v, int max_iter, float epsilon, float alfa) // v are directions, should be normalized
    {

        float[][] v_conv = new float[v.size()][2];
        for (int i = 0; i < v.size(); i++) {
            v_conv[i][0] = v.get(i)[0]; // should be normalized
            v_conv[i][1] = v.get(i)[1];
        }

        // auxiliary variable for iteration (to avoid allocation inside the loop)
        float[] new_v = new float[2];

        for (int i = 0; i < v_conv.length; i++) {
            int iter = 0;
            double d;

            do {

                runOne(v_conv[i], new_v, v, alfa); // new_v is normalized

                float dot_prod = new_v[0] * v_conv[i][0] + new_v[1] * v_conv[i][1];
                dot_prod = (dot_prod>1)? 1 : dot_prod;
                d = Math.acos(dot_prod);

                v_conv[i][0] = new_v[0];
                v_conv[i][1] = new_v[1];

                iter++;
            }
            while (iter < max_iter&& d > epsilon);

        }
        return v_conv;
    }

    private static void runOne(float[] curr_v, float[] new_v, ArrayList<float[]> v, float alfa)
    {

        float sum 	= 0;
        new_v[0] 	= 0;
        new_v[1] 	= 0;

        for (int l = 0; l < v.size(); l++) { // loop all directions
            if (Math.acos(curr_v[0]*v.get(l)[0] + curr_v[1]*v.get(l)[1]) <= alfa) { // if they are unit vectors it is not necessary to divide with norms here
                sum += 1;
                new_v[0] += 1 * v.get(l)[0];
                new_v[1] += 1 * v.get(l)[1];
            }
        }

        new_v[0] /= sum;
        new_v[1] /= sum;
        float norm = (float) Math.sqrt(Math.pow(new_v[0],2)+Math.pow(new_v[1],2));

        if (sum>0 && norm>0) {
            // normalize (because we're mean-shifting directions)
            new_v[0] /= norm;
            new_v[1] /= norm;
        }
        else {
            new_v[0] = curr_v[0];
            new_v[1] = curr_v[1];
        }

    }

    private static int[] clustering(float[][] values, float threshold_dists)
    {
        // indxs represent indexes of values that need to be clustered according to their values read before
        // dists are the distances (idxs.length * idxs.length),
        // threshold_dists is the distance limit
        // output is list of unique labels

        int[] labels = new int[values.length];
        for (int i = 0; i < labels.length; i++) labels[i] = i; // each object gets its unique label

//		System.out.println("BEFORE:");
//		System.out.println(Arrays.toString(labels));

        // not really efficient but calculate it here
        float[][] dists = new float[values.length][values.length];
        for (int i = 0; i < values.length; i++) {
            for (int j = i; j < values.length; j++) {
                if(i==j) {
                    dists[i][j] = 0;
                }
                else {

                    float dot_prod = values[i][0]*values[j][0] + values[i][1]*values[j][1]; // vi * vj
                    dot_prod = (dot_prod>1)? 1 : dot_prod;
                    float dij = (float) Math.acos(dot_prod);
                    dists[i][j] = dij;
                    dists[j][i] = dij;
                }
            }
        }
//		for (int i =0 ; i<dists.length; i++) IJ.log(""+Arrays.toString(dists[i]));

        for (int i = 0; i < values.length; i++) {

            // one versus the rest
            for (int j = 0; j < values.length; j++) {

                if (i != j) {

                    if (dists[i][j]<=threshold_dists) {

                        if (labels[j] != labels[i]) {

                            int currLabel = labels[j];
                            int newLabel  = labels[i];

                            labels[j] = newLabel;

                            //set all that also were currLabel to newLabel
                            for (int k = 0; k < labels.length; k++)
                                if (labels[k]==currLabel)
                                    labels[k] = newLabel;

                        }

                    }

                }

            }

        }

//		System.out.println("AFTER:");
//		System.out.println(Arrays.toString(labels));

        return labels;

    }

    private static ArrayList<float[]> extracting(int[] labels, ArrayList<float[]> vals, int min_count)
    {

        boolean[] checked = new boolean[labels.length];

        ArrayList<float[]> out = new ArrayList<float[]>();
        ArrayList<Integer> cnt = new ArrayList<Integer>();    		// to make sure that it outputs sorted list

        for (int i = 0; i < labels.length; i++) {
            if (!checked[i]) {

                float centroid_x = vals.get(i)[0]; // idxs[i]
                float centroid_y = vals.get(i)[1]; // idxs[i]
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            centroid_x += vals.get(j)[0]; // idxs[j]
                            centroid_y += vals.get(j)[1]; // idxs[j]
                            count++;
                            checked[j] = true;

                        }
                    }
                }

                if (count >= min_count) {
                    out.add(new float[]{centroid_x/count, centroid_y/count});
                    cnt.add(count);
                }

            }
        }


        // print before sorting
//		for (int ii = 0; ii < cnt.size(); ii++) {
//			IJ.log(ii + " : " + cnt.get(ii) + " points,  at " + Arrays.toString(out.get(ii)));
//		}

        // sort by the counts (take from Sorting.java) descending indexes
        int[] desc_idx = Tools.descending(cnt);      // it will change cnt list as well!!!!
        int clusters_nr = (desc_idx.length>4)? 4 : desc_idx.length ;

        ArrayList<float[]> out_sorted = new ArrayList<float[]>(clusters_nr); // top 4  if there are as many

//		out_sorted = new float[clusters_nr][2];
//		System.out.println(clusters_nr+" clusters allocated");
//		int[] clustered_counts = new int[clusters_nr];

        for (int ii=0; ii<clusters_nr; ii++) {

//			_clustered_directions[ii][0] = out.get(desc_idx[ii])[0];
//			_clustered_directions[ii][1] = out.get(desc_idx[ii])[1];
//			clustered_counts[ii] = cnt.get(ii);

            float vx 	= out.get(desc_idx[ii])[0];
            float vy 	= out.get(desc_idx[ii])[1];
            float vcnt 	= cnt.get(ii);  // because cnt is already sorted

            out_sorted.add(new float[]{vx, vy, vcnt}); // add top 1,2,3 or 4 directions based on the count
        }

        return out_sorted;

    }

    private static float Deg2Rad(float ang_deg)
    {return (float) ((ang_deg/180f)*Math.PI);}
}
