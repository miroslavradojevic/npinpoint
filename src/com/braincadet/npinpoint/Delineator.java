package com.braincadet.npinpoint;

import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.process.ByteProcessor;

import java.awt.*;
import java.util.ArrayList;
import java.util.Arrays;

/**  */
public class Delineator extends Thread {

    // associate the peaks and link follow-up points, delineate local structure, extract smoothness feature

    private int begN, endN;

    // INPUT
    public static int[][] 	    i2xy;                       // selected locations
    public static int[][]     	xy2i;                       // need for recursion
    public static int[][]       peaks_i;             	    // list of extracted peaks: N x 4 every FG location with 4 extracted peaks in indexed format
    public static int[][]		peaks_w;                    // weight assigned to each peak (for expansion)
    public static float[][]		inimg_xy;				    // input image (necessary for feature extraction)
    public static boolean[] 	critpoint_candidate;        // if profile had enough of range (same length as i2xy)

    public static float     D;								// neuron diameter parameter
    public static int       M;                              // how much it expands recursively from the center
    public static float     minCos;                         // allowed derail
    public static int		L = Integer.MIN_VALUE;          // will define how many are taken along the diameter, in radial direction

    private static int 		windowSize = 5;                 // vxy (refined vectors) calculate by averaging neighbouring directions
    public static float 	samplingStep    = .7f;          // cross-side sampling step
    private static float    samplingStepLongitudinal = Integer.MIN_VALUE; // sampling along the streamline of patches

    public static int       dim;                            // cross profile length (number of samples in a cross-profile)
    public static int      	dim_half;                       // cross profile half-length
    public static int 		half_window;

    private static float	clusteringDisk = Float.NaN;


    private static int      Nimg = Integer.MIN_VALUE;// wrt sphere size
    private static float    outer_radius = Float.NaN;

    // OUTPUT
    public static int[][][]     delin2;                     // N(foreground locs.) x 4(max. threads) x (1..M) (follow-up locs.) contains index for each location
    public static float[][][][] xy2;        				// N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
    public static boolean[][][] xy2sel;						// N                   x 4                   x ((1..M) x L)
    public static float[][][][] vxy2;        				// N(foreground locs.) x 4(max. threads) x 2 x ((1..M) x L) (2 dimensional)
    public static float[][]		smoothness;					// N(foreground locs.) x 4(max. threads) smoothness score

    public Delineator(int n0, int n1)
    {
        this.begN = n0;
        this.endN = n1;
    }

    public static void loadTemplate(
            int[][] 		_i2xy,
            int[][] 		_xy2i,
            boolean[] 	_critpoint_candidate,
            int[][] 		_peaks_i,
            int[][] 		_peaks_w,
            float[][] 	_inimg_xy,
            Sphere2D		_extractionSphere,
            int     		_M,
            float   		_minCos
    )
    {

        i2xy = _i2xy;
        xy2i = _xy2i;
        peaks_i = _peaks_i;
        peaks_w = _peaks_w;
        inimg_xy = _inimg_xy;
        critpoint_candidate = _critpoint_candidate;

        D 				= _extractionSphere.getNeuronDiameter();
        M               = _M;
        minCos          = _minCos;
        L               = (int) (Math.ceil(_extractionSphere.getNeuronDiameter()) + 3);  // added two because of second derivative calcualtion on the refined locations

        samplingStepLongitudinal = D / (float)(L-1);

//		System.out.println("longitudinal: " + samplingStepLongitudinal);

        // dim for the cross-profile
        dim_half = (int) Math.ceil( D / (samplingStep*2) );
        dim = 2*dim_half + 1;
        half_window = (int) Math.round(((float)dim/4f - 1) * 0.5); // dim for the cross-profile window
//		System.out.println("" + dim + " ,,,, " + half_window);
        half_window = (half_window<1)? 1 : half_window;
//		System.out.println("w2 = " + half_window);

        clusteringDisk = 1.5f * samplingStepLongitudinal;

        outer_radius = _extractionSphere.getOuterRadius();
        Nimg = (int) Math.ceil(outer_radius);

        // allocate output
        delin2 				= new int[i2xy.length][4][M];
        // initialize with -2 (same meaning as with peaks_i)
        for (int i=0; i<delin2.length; i++)
            for (int j=0; j<delin2[i].length; j++)
                for (int k=0; k<delin2[i][j].length; k++)
                    delin2[i][j][k] = -2;

        xy2 		= new float[i2xy.length][4][2][];      // N (foreground points) x 4 (branches) x 2 (x,y) x M*1..L (along branches) (x,y)
        xy2sel      = new boolean[i2xy.length][4][];       // N (foreground points) x 4 (branches) x M*1..L
        smoothness 	= new float[xy2.length][4];            // N(foreground locs.)   x 4 (branches)
        vxy2 		= new float[i2xy.length][4][2][];      // N (foreground points) x 4 (branches) x 2 (vx,vy) x M*1..L (along branches) (vx,vy)

    }

    public void run()
    {

        //*** auxiliary
        // to store selected local cross section xy values (side output when calculating refined locs)
        float[][][] xy_local = new float[4][2][]; // same variable will be used for all the locations

        ArrayList<ArrayList<float[]>> biggest_clusters = new ArrayList<ArrayList<float[]>>(4); // allocate 4 clusters
        for (int i = 0; i < 4; i++) biggest_clusters.add(new ArrayList<float[]>(1));  // size = 4x0

        byte[][][] mask = new byte[4][2*Nimg+1][2*Nimg+1];

        ArrayList<ArrayList<float[]>> selected_xy_local = new ArrayList<ArrayList<float[]>>(4);
        for (int i = 0; i < 4; i++) selected_xy_local.add(new ArrayList<float[]>(1)); // size = 4X0

        // cross-profile allocation here
        float[]  	patch_profile 	= new float[dim]; // patch_profile used for cross-sections - initialize, auxiliary variable filled up at each cross-section when calculating refined xy points
        float[]		cumsum 			= new float[dim]; // cummulative sum (for local average calculation)


        //*** auxiliary

        for (int locationIdx = begN; locationIdx < endN; locationIdx++) {

            if (!critpoint_candidate[locationIdx]) { // ruled out because of insignificant profile variation

                delin2[locationIdx] = null;
                xy2[locationIdx] = null;
                vxy2[locationIdx] = null;

            }
            else {

                /**************************************************************************
                 delin2[locationIdx]
                 **************************************************************************/
                for (int pp = 0; pp<peaks_i[locationIdx].length; pp++) {  // loop 4 allocated branches (1st generation)
                    // access individual peaks at this point (-2:not exist, -1:background, >=0:foreground)
                    if (peaks_i[locationIdx][pp] >= 0) {

                        //if the peak exists in the foreground
                        delin2[locationIdx][pp][0] = peaks_i[locationIdx][pp];

                        int curr_index, prev_index, next_index;
                        curr_index = peaks_i[locationIdx][pp];
                        prev_index = locationIdx;

                        for (int m=1; m<M; m++) { // follow the recursion for the rest of the indexes

                            next_index = getNext(prev_index, curr_index); // recursion : prev+curr->next index
                            // >=0 if there was at least one in foreground
                            // will give -1 if the next one was in the background and there was no other foreground peak
                            // -2 if it there was no other foreground peak and background either

                            if (next_index>=0) {
                                delin2[locationIdx][pp][m] = next_index;     // store it in output matrix
                            }
                            else if (next_index==-1){ // not found but in the background
                                delin2[locationIdx][pp][m] = -1;
                                // fill the rest with -1 once one fell out
                                for (int m_aux=m+1; m_aux<M; m_aux++) delin2[locationIdx][pp][m_aux] = -1;
                                break; // stop looping further along the streamline here
                            }
                            else if (next_index==-2) { // not found at all
                                delin2[locationIdx][pp][m] = -2;
                                // fill the rest with -2 once one had a dead end
                                for (int m_aux=m+1; m_aux<M; m_aux++) delin2[locationIdx][pp][m_aux] = -2;
                                break; // stop looping further along the streamline here
                            }

                            // recursion - to expand further in the next iteration
                            prev_index = curr_index;
                            curr_index = next_index;

                        }

                    }
                    else if (peaks_i[locationIdx][pp] == -1) {

                        // 1st generation peak index was -1 (peak in the background according to the mask)
                        // this is a case of an incomplete delineation, having it's parts in the background
                        // it can be anything and mask should not influence the final decision making
                        // if the peak was discarded as it did not exist - wrong it can affect the decision
                        // solution: keep only those delineations unaffected by the mask, where the whole delineation is in foreground
                        // assign the whole delineation as incomplete here - giving -1 to all elements and finish
                        for (int m=0; m<M; m++) delin2[locationIdx][pp][m] = -1;

                    }
                    else if (peaks_i[locationIdx][pp] == -2) {

                        // those are the places that were not filled up with peaks
                        // propagate it to delin2
                        for (int m=0; m<M; m++) delin2[locationIdx][pp][m] = -2;

                    }

                }

                /**************************************************************************
                 examples delin2[locationIdx][4][M], (M=3):
                 23  34 56
                 678 -1 -1      -> went into background or couldn't find successor
                 -2  -2 -2      -> does not exist
                 **************************************************************************/

                /**************************************************************************
                 xy2[locationIdx]
                 xy_local[4][2][]  used to accumulate patch variables (xy) for smoothness calculation
                 **************************************************************************/
                for (int b = 0; b<delin2[locationIdx].length; b++) { // loop 4 branches

                    if (delin2[locationIdx][b][0]==-1) {
                        // whole streamline is missing: the first one was -1, recursion was stopped
                        xy2[locationIdx] = null;
                        break; 	// stop looping branches further
                    }
                    else if (delin2[locationIdx][b][0]==-2) {
                        // no streamline here
                        xy2[locationIdx][b] = null;
                        xy_local[b] = null;
                    }
                    else if (delin2[locationIdx][b][0]>=0) {
                        // there is at least one patch in the streamline
                        // loop the rest to count how many patches there are to allocate the array
                        int count_patches = 1;
                        for (int m=1; m<M; m++) {
                            if (delin2[locationIdx][b][m] >= 0) count_patches++;
                            else break; // because the rest are filled up with -1 or -2 anyway
                        }

                        // allocate (now we know how much to allocate)
                        xy2[locationIdx][b][0] 	= new float[count_patches*L]; // x coordinates allocate
                        xy2[locationIdx][b][1] 	= new float[count_patches*L]; // y coordinates allocate
                        xy_local[b] = new float[2][count_patches*L] ;   // auxiliary variable allocate

                        // fill  up the allocated arrays for each patch
                        for (int m = 0; m<M; m++) {      					// loop patches outwards, from the beginning

                            if (delin2[locationIdx][b][m]>=0) {  // will be at least one true because it was used to count the patches in previous loop

                                // there is a patch, add the features to the matrix
                                int curr_i = delin2[locationIdx][b][m];
                                int curr_x = i2xy[curr_i][0];
                                int curr_y = i2xy[curr_i][1];

                                int prev_i, prev_x, prev_y;

                                if (m==0) {
                                    prev_x = i2xy[locationIdx][0];  // central location
                                    prev_y = i2xy[locationIdx][1];
                                }
                                else{
                                    prev_i = delin2[locationIdx][b][m-1];
                                    prev_x = i2xy[prev_i][0];
                                    prev_y = i2xy[prev_i][1];
                                }

                                // get refined locations sampled from the local patch (aligned with the patch)
                                localPatchRefinedLocs(
                                        prev_x, prev_y, curr_x, curr_y,
                                        m * L,        // start from
                                        xy2[locationIdx][b],
                                        xy_local[b],
                                        patch_profile,  // one profile is allocated and referenced here, filled up each time
                                        cumsum
                                ); //side output: local values for derivation are stored in xy_local[4][2][]
                            }
                            else break; // because the rest are filled up with -1 or -2 anyway

                        }  // xy2 is formed for this branch

                    } // branch case

                } // loop branches

                // takes into account the spatioal distribution of the refined points along branches
                /**************************************************************************
                 biggest_clusters (local variable) using xy2
                 **************************************************************************/
                extractClustersForRefinedStreamlineLocations(biggest_clusters, xy2[locationIdx]); // expects refined points in array [locId][4][2][M*(1..L)]

                /**************************************************************************
                 mask using biggest_clusters
                 **************************************************************************/
                extractLocalClusterMask(mask, biggest_clusters, locationIdx); // forms mask that marks spatial regions corresponding to the biggest cluster for each branch

                /**************************************************************************
                 xy2sel[locationIdx]
                 selected_xy_local (list judging on whether the corresponding xy2 belonged to the biggest cluster in other branch)
                 **************************************************************************/
                extractRefinedStreamlineSelection(xy2sel[locationIdx], selected_xy_local, xy2[locationIdx], locationIdx, xy_local, biggest_clusters, mask); // extract selection to calculate second derivative (hence smoothness)

                /**************************************************************************
                 smoothness[locationIdx] using selected_xy_local
                 **************************************************************************/

                for (int i = 0; i < selected_xy_local.size(); i++) {

                    int nr_pts = selected_xy_local.get(i).size();

                    if (nr_pts>=(float)L/2) { // 3 are enough for second derivative calculation, needs more for robustness

                        smoothness[locationIdx][i] = 0;
                        for (int j = 0; j <nr_pts-2; j++) {

                            float x0 = selected_xy_local.get(i).get(j)[0];
                            float y0 = selected_xy_local.get(i).get(j)[1];

                            float x1 = selected_xy_local.get(i).get(j+1)[0];
                            float y1 = selected_xy_local.get(i).get(j+1)[1];

                            float x2 = selected_xy_local.get(i).get(j+2)[0];
                            float y2 = selected_xy_local.get(i).get(j+2)[1];

                            float f2 = ( (y2-y1)/(x2-x1) - (y1-y0)/(x1-x0) ) / (x1-x0);
                            float dx = x1-x0;

                            smoothness[locationIdx][i] += f2*f2*dx;

                        }

                    }
                    else {
                        smoothness[locationIdx][i] = Float.NaN;
                    }

                }

                /**************************************************************************
                 vxy2[locationIdx]
                 **************************************************************************/
                if (xy2[locationIdx]!=null) {
                    for (int b = 0; b<xy2[locationIdx].length; b++) {
                        if (xy2[locationIdx][b]!=null) {

                            int to_allocate = xy2[locationIdx][b][0].length;
                            vxy2[locationIdx][b][0] = new float[to_allocate]; // vx
                            vxy2[locationIdx][b][1] = new float[to_allocate]; // vy

                            for (int l=0; l<to_allocate; l++) {

                                float avg_v_x=0, avg_v_y=0;

                                for (int l_nbr=l-windowSize/2; l_nbr<=l+windowSize/2; l_nbr++) {

                                    if (l_nbr>=0 && l_nbr+1<to_allocate) {

                                        float x1 = xy2[locationIdx][b][0][l_nbr];
                                        float y1 = xy2[locationIdx][b][1][l_nbr];
                                        float x2 = xy2[locationIdx][b][0][l_nbr+1];
                                        float y2 = xy2[locationIdx][b][1][l_nbr+1];

                                        float norm = (float) Math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
                                        float vx = (x2-x1)/norm;
                                        float vy = (y2-y1)/norm;

                                        avg_v_x += vy;
                                        avg_v_y += -vx;

                                    }

                                }

                                float norm = (float) Math.sqrt(avg_v_x*avg_v_x + avg_v_y*avg_v_y);
                                if (norm>Float.MIN_VALUE) {
                                    vxy2[locationIdx][b][0][l] = avg_v_x/norm;
                                    vxy2[locationIdx][b][1][l] = avg_v_y/norm;
                                }
                                else {
                                    vxy2[locationIdx][b][0][l] = 0;
                                    vxy2[locationIdx][b][1][l] = 0;
                                }

                            }

                        }
                        else {
                            vxy2[locationIdx][b] = null;
                        }
                    }
                }
                else {
                    vxy2[locationIdx] = null;
                }

            } // it was critpoint candidate (selected from foreground points)

        }

    }

    private static int getNext(int prev_index, int curr_index)
    {

        // these are stacked as XY - take care on that
        int prevX = i2xy[prev_index][0];    // X
        int prevY = i2xy[prev_index][1];    // Y

        int currX = i2xy[curr_index][0];    // X
        int currY = i2xy[curr_index][1];    // Y

        // check peaks at curr
        int[]   pks4xI  = peaks_i[curr_index]; // list all indexes
        int[] 	pks4xW	= peaks_w[curr_index];

        // take the one with highest weight that is above minCos
        int next_index 	= -2;
        int next_weight = Integer.MIN_VALUE;
        boolean found_valid_followup = false;

        for (int p = 0; p<pks4xI.length; p++) {

            if (pks4xI[p]>=0) {
                // peak in foreground
                int check_I = pks4xI[p];
                int next_X 	= i2xy[check_I][0];
                int next_Y	= i2xy[check_I][1];

                double cosAng =
                        (
                                (currX-prevX)*(next_X-currX) + (currY-prevY)*(next_Y-currY)           // + (currZ-prevZ)*(nextZ-currZ)
                        )
                                /
                                (
                                        Math.sqrt( Math.pow(currX-prevX, 2) + Math.pow(currY-prevY, 2) ) *  //  + Math.pow(currZ-prevZ, 2)
                                                Math.sqrt( Math.pow(next_X-currX, 2) + Math.pow(next_Y-currY, 2) )    //  + Math.pow(nextZ-currZ, 2)
                                );

                if (cosAng>minCos) {
                    // peak is outwards pointing
                    if (pks4xW[p]>next_weight) {
                        found_valid_followup = true;
                        next_weight = pks4xW[p];
                        next_index 	= pks4xI[p];
                    }
                }

            }
            else {
                // peak in the background and there is nothing acceptable in foreground so far
                if (!found_valid_followup) {
                    next_index = -1;
                    //next_weight neutral
                }
                else {
                    // there is a follow up to use, ignore background point
                }

            }
        }

        return next_index;
        // will give -1 if the next one was in the background and there was no other foreground peak
        // -2 if it there was no other foreground peak and background either
        // >=0 if there was at least one in foreground

    }

    private static void localPatchRefinedLocs(float x1, float y1, float x2, float y2,
                                              int init_index,
                                              float[][] refined_centerline_locs_xy,
                                              float[][] refined_local_ptch_xy,
                                              float[] _patch_profile,
                                              float[] _cumsum
    )
    {
        /*
            standard way to loop through a patch defined with 2 points
         */
        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        boolean found_local_max;
        float cross_profile_max;
        float loc_avg;

        for (int ii=0; ii<L; ii++) {  // loops L of them in radial direction with 0(root coordinate) and L included

            // patch_profile fill-up
            for (int jj=-dim_half; jj<=dim_half; jj++) {
                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;
                _patch_profile[jj+dim_half] = Interpolator.interpolateAt(curr_x, curr_y, inimg_xy);
                _cumsum[jj+dim_half] = _patch_profile[jj+dim_half];
                if (jj+dim_half>0)
                    _cumsum[jj+dim_half] += _cumsum[jj+dim_half-1];
            }
            // patch_profile, cumsum filled-up

            // find max, max is calculated as local average
            found_local_max = false;
            cross_profile_max = Float.NEGATIVE_INFINITY;

            // ii marks longitudinal index, tangential index is set to
//			float refined_x = x_root + ii * samplingStepLongitudinal * vx + (  prev_local_max_idx  ) * samplingStep * vy;
//			float refined_y = y_root + ii * samplingStepLongitudinal * wx + (  prev_local_max_idx  ) * samplingStep * wy;
//			float yptch = (  prev_local_max_idx  ) * samplingStep; // starts from the patch border, keep patch cross-projection along the cross-section

//			offset_min = Float.MAX_VALUE;//Math.abs(1 - (patch_profile.length-1)/2) * samplingStep;
//            int index_from_center   = -1;
//            int index_from_beg      = -1;

            float refined_x = Float.NaN, refined_y=Float.NaN, yptch=Float.NaN;  // there will be stored

            for (int jj=0; jj<dim; jj++) {

                if (jj<=half_window)
                    loc_avg = _cumsum[jj + half_window] / (1 + half_window + jj);
                else if (jj>half_window && jj<=dim-1-half_window)
                    loc_avg = (_cumsum[jj + half_window] - _cumsum[jj - 1 - half_window]) / (1 + 2 * half_window);
                else
                    loc_avg = (_cumsum[dim - 1] - _cumsum[jj - 1 - half_window]) / (1 + half_window + (dim - 1 - jj));

                if (loc_avg>cross_profile_max) {

                    cross_profile_max = loc_avg;

                    int index_distance_from_center = jj-dim_half;

                    refined_x  = x_root + ii * samplingStepLongitudinal * vx + index_distance_from_center * samplingStep * vy;
                    refined_y  = y_root + ii * samplingStepLongitudinal * wx + index_distance_from_center * samplingStep * wy;
                    yptch      = index_distance_from_center * samplingStep;

                    found_local_max = true;

                }

            }

//			if (!found_local_max) System.out.println("couldn't optimize cross-section");

            // store refined_x, refined_y, yptch
            refined_centerline_locs_xy[0][init_index + ii] = refined_x; // ii is radial index
            refined_centerline_locs_xy[1][init_index + ii] = refined_y;
            // store cross y values
            refined_local_ptch_xy[0][init_index + ii] = ii * samplingStepLongitudinal;
            refined_local_ptch_xy[1][init_index + ii] = yptch;

        }  // end looping radial indexes

    }

    public static ImageStack getSmoothnessDistribution(int nr_bins)
    {

        ArrayList<Float> concatenate_smoothness_scores = new ArrayList<Float>(smoothness.length*4);

        System.out.println(smoothness.length + " is amount of points");

        int cnt = 0;

        for (int i = 0; i < smoothness.length; i++) {

            if (smoothness[i]!=null) {

                for (int j = 0; j < smoothness[i].length; j++) {

                    if (smoothness[i][j]!=Float.NaN) {
                        concatenate_smoothness_scores.add(smoothness[i][j]);

                    }
                    else {
                        cnt++;
                    }


                }

            }


        }

        System.out.println("total " + concatenate_smoothness_scores.size() + " smoothness scores");
        System.out.println("minimum = " + Hist.find_min(concatenate_smoothness_scores));
        System.out.println("maximum = " + Hist.find_max(concatenate_smoothness_scores));

        ImageStack is_out = new ImageStack(528, 255);

        float[] bins = new float[nr_bins];
        float[] distribution = new float[nr_bins];

        Hist.getDistribution(concatenate_smoothness_scores, nr_bins, bins, distribution);
        Plot p = new Plot("", "", "", bins, distribution, Plot.LINE);
        is_out.addSlice("", p.getProcessor());
        return is_out;

    }

    public static ImageStack getClusterMask(int atX, int atY)
    {

        int locIdx = xy2i[atX][atY];
        ImageStack is_out = new ImageStack(2*Nimg+1, 2*Nimg+1);


        if (locIdx>=0) {

            // define clusters
            ArrayList<ArrayList<float[]>> biggest_clusters = new ArrayList<ArrayList<float[]>>(4); // allocate 4 clusters
            for (int i = 0; i < 4; i++) biggest_clusters.add(new ArrayList<float[]>(1));  // size = 4x0

            // extract clusters
            extractClustersForRefinedStreamlineLocations(biggest_clusters, xy2[locIdx]);

            // define mask
            byte[][][] mask = new byte[4][2*Nimg+1][2*Nimg+1];

            // extract mask
            extractLocalClusterMask(mask, biggest_clusters, locIdx);

            // reformat the array
            int W = mask[0].length;
            int H = mask[0][0].length;
            int L = mask.length;

            byte[][] mask_out = new byte[L][W*H];
            for (int cmask = 0; cmask < L; cmask++) {
                for (int xmask = 0; xmask < W; xmask++) {
                    for (int ymask = 0; ymask < H; ymask++) {
                        mask_out[cmask][ymask*W+xmask] = mask[cmask][xmask][ymask];
                    }
                }
            }

            // export current mask
            for (int i = 0; i < mask_out.length; i++) {
                ByteProcessor bp_out = new ByteProcessor(W, H, mask_out[i]);
                is_out.addSlice(bp_out);
            }


//			System.out.println("Extracted biggest clusters at ("+atX+","+atY+"): ");
//			for (int i = 0; i < biggest_clusters.size(); i++) {
//				for (int j = 0; j < biggest_clusters.get(i).size(); j++) {
//					System.out.println("cluster " + i + " ,element " + j + " : " + Arrays.toString(biggest_clusters.get(i).get(j)));
//				}
//			}
//
//			System.out.println("smoothness : "+Arrays.toString(smoothness[locIdx]));

        }

        return is_out;

    }

    public static Overlay getDelineationOverlay(int atX, int atY)
    {

        // return the delineated local structure Overlay for visualization
        Overlay ov = new Overlay();

        /*
            adds generations of peaks in green
         */
        PeakExtractor.getPeaks(atX, atY, 3, ov);   // peaks will be stored in ov

        float Rd = 1.5f; // radius of the circles written for delineation
        Color cd = Color.RED;
        float wd = 0.25f;

		/*
		 central location
		  */
        OvalRoi ovalroi = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
        ovalroi.setFillColor(cd);
        ovalroi.setStrokeWidth(wd);
        ov.add(ovalroi);

        /*
            add local cross-section locations (re-calculate) and patch delineation spots (delin2)
         */
        if (xy2i[atX][atY]>=0) {

            int[][] delin_at_loc = delin2[xy2i[atX][atY]];

            if (delin_at_loc!=null) {

                for (int b = 0; b<delin_at_loc.length; b++) {           // loop 4 branches

                    for (int m=0; m<M; m++) {

                        if (delin_at_loc[b][m] >= 0) {

                            int pt_idx = delin_at_loc[b][m]; // there is a point to add
                            int pt_x = i2xy[pt_idx][0];
                            int pt_y = i2xy[pt_idx][1];

                            ovalroi = new OvalRoi(pt_x-(Rd/2)+.5f, pt_y-(Rd/2)+.5f, Rd, Rd); // add the point to the overlay
                            ovalroi.setStrokeColor(cd);
                            ovalroi.setStrokeWidth(wd);
                            ov.add(ovalroi);

                            // find previous indexes
                            int prev_i, prev_x, prev_y;

                            if (m==0) {
                                prev_x = atX;
                                prev_y = atY;
                            }
                            else{
                                prev_i = delin_at_loc[b][m-1];
                                prev_x = i2xy[prev_i][0];
                                prev_y = i2xy[prev_i][1];
                            }

                            ArrayList<OvalRoi> line_pts = localPatchCrossProfilesLocs(prev_x, prev_y, pt_x, pt_y);  // will be recalculated
                            for (int aa = 0; aa < line_pts.size(); aa++) ov.add(line_pts.get(aa));

                        }

                    }

                }

            } // delin2 != null
            else{
                // do nothing delineatoin was null here because profile was too flat
            }

        }

        /*
        add refined streamline locations (read from xy2)
         */
        int locationIdx = xy2i[atX][atY];

        if (locationIdx!=-1) {

            if (xy2[locationIdx]!=null) {

                for (int b = 0; b<xy2[locationIdx].length; b++) {

                    // check whether there is a refinement at all here, add in case there is
                    if (xy2[locationIdx][b]!=null) {


                        // loop points to add them
                        if (xy2[locationIdx][b][0]!=null) {

                            for (int l=0; l<xy2[locationIdx][b][0].length; l++) {

                                float refined_x = xy2[locationIdx][b][0][l];
                                float refined_y = xy2[locationIdx][b][1][l];

                                ovalroi = new OvalRoi(refined_x-(samplingStep/2)+.5f, refined_y-(samplingStep/2)+.5f, samplingStep, samplingStep);
                                ovalroi.setFillColor(xy2sel[locationIdx][b][l]?Color.YELLOW:new Color(1f,1f,0f,0.5f));     // Color.YELLOW
                                ovalroi.setStrokeWidth(samplingStep/4);
                                ov.add(ovalroi);

                            }


                        }
                        else {

                            System.out.println("xy2 null");


                            float old_Rd = Rd;

                            Rd = 3*Rd;
                            OvalRoi ovalroi_test = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
                            ovalroi_test.setFillColor(Color.MAGENTA);
                            ovalroi_test.setStrokeWidth(3);
                            ov.add(ovalroi_test);
                            System.out.println("ADDED");


                            Rd = old_Rd;

                        }



                    }

                }
            }
            else {
                // do nothing, do nothing, points were not extracted because profile was flat
            }

        }

        /*
        add lines marking the refined cross sections
         */
        if (locationIdx!=-1) {

            if (vxy2[locationIdx]!=null) {

                for (int b=0; b<vxy2[locationIdx].length; b++) {

                    if (vxy2[locationIdx][b]!=null) {






                        // loop points to add them
                        if (vxy2[locationIdx][b][0]==null) {
                            System.out.println("vxy2 null");

                            float old_Rd = Rd;

                            Rd = 3*Rd;
                            OvalRoi ovalroi_test = new OvalRoi(atX-(Rd/2)+.5f, atY-(Rd/2)+.5f, Rd, Rd);
                            ovalroi_test.setStrokeColor(Color.MAGENTA);
                            ovalroi_test.setStrokeWidth(3);
                            ov.add(ovalroi_test);

                            Rd = old_Rd;

                        }
                        else {

                            for (int l=0; l<vxy2[locationIdx][b][0].length; l++) {

                                float end_1_x = xy2[locationIdx][b][0][l] - dim_half * samplingStep * vxy2[locationIdx][b][0][l];
                                float end_1_y = xy2[locationIdx][b][1][l] - dim_half * samplingStep * vxy2[locationIdx][b][1][l];

                                float end_2_x = xy2[locationIdx][b][0][l] + dim_half * samplingStep * vxy2[locationIdx][b][0][l];
                                float end_2_y = xy2[locationIdx][b][1][l] + dim_half * samplingStep * vxy2[locationIdx][b][1][l];

                                Line lne = new Line(end_1_x+.5f, end_1_y+.5f, end_2_x+.5f, end_2_y+.5f);
                                lne.setStrokeColor(Color.CYAN);
                                ov.add(lne);

                            }


                        }












                    }

                }

            }
            else {
                // profile was too flat here to extract
            }

        }

        return ov;

    }

    private static ArrayList<OvalRoi> localPatchCrossProfilesLocs(float x1, float y1, float x2, float y2)
    {

        float l = (float) Math.sqrt(Math.pow(x2-x1, 2)+Math.pow(y2-y1, 2));
        float vx = (x2-x1)/l;
        float vy = (y2-y1)/l;
        float wx = vy;
        float wy = -vx;

        float x_root = x2 - vx * D;
        float y_root = y2 - wx * D;

        float R = samplingStep/2;

        ArrayList<OvalRoi> pts = new ArrayList<OvalRoi>(dim*L);

        for (int ii=0; ii<L; ii++) { // loops L of them in radial direction

            for (int jj=-dim_half; jj<=dim_half; jj++) { // loops vector w

                float curr_x = x_root + ii * samplingStepLongitudinal * vx + jj * samplingStep * vy;
                float curr_y = y_root + ii * samplingStepLongitudinal * wx + jj * samplingStep * wy;

                OvalRoi pt = new OvalRoi(curr_x+.5f-(R/2), curr_y+.5f-(R/2), R, R);
                pt.setFillColor(Color.BLUE);
                pt.setStrokeWidth(R/2);
                pts.add(pt);

            }

        }

        return pts;

    }

    private static ArrayList<float[]> largest_cluster(float[] disks_x, float[] disks_y, float disk_r) // [x, y, r] all the disks will have the same radius
    {

        int[] labels = new int[disks_x.length];
        for (int i = 0; i < labels.length; i++) labels[i] = i;

        for (int i = 0; i < disks_x.length; i++) {

            // one versus the rest
            for (int j = 0; j < disks_x.length; j++) {

                if (i!=j) {

                    double dst2 	= Math.pow(disks_x[i]-disks_x[j], 2) + Math.pow(disks_y[i]-disks_y[j], 2);
                    double rd2 		= Math.pow(disk_r+disk_r, 2);

                    if (dst2<=rd2) {  // they are neighbours

                        if (labels[j]!=labels[i]) {

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

        boolean[] checked = new boolean[labels.length];

        // output list of list of locations clustered together
        ArrayList<float[]> lagrest_cluster = new ArrayList<float[]>();
        ArrayList<float[]> cluster_xy      = new ArrayList<float[]>();

        int max_size = Integer.MIN_VALUE;

        for (int i = 0; i < labels.length; i++) {

            if (!checked[i]) {

                cluster_xy.clear(); // reset with each new label
                cluster_xy.add(new float[]{disks_x[i], disks_y[i]});
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < labels.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            cluster_xy.add(new float[]{disks_x[j], disks_y[j]});
                            checked[j] = true;

                        }
                    }
                }

                if (cluster_xy.size()>max_size) {
                    // set it as output cluster (maybe there is a more efficient way to do this)
                    lagrest_cluster.clear();
                    for (int cc = 0; cc < cluster_xy.size(); cc++) {
                        lagrest_cluster.add(cluster_xy.get(cc).clone());
                    }
                }
//                clusters.add(cluster_xy);

            }
        }

        // take the largest one

        return lagrest_cluster; // cluster labels for each disc

    }

    private static boolean belongs_to_cluster(float point_x, float point_y, ArrayList<float[]> cluster_xy, float disk_r)
    {
        // loop the contents of the cluster & see if the point is close to any of the cluster points
        // if it is - say it belongs
        boolean it_belongs = false;
        for (int i = 0; i < cluster_xy.size(); i++) {
            if (dist2(point_x, point_y, cluster_xy.get(i)[0], cluster_xy.get(i)[1]) <= disk_r*disk_r) {
                it_belongs = true;
                break; // one is enough to know
            }
        }
        return it_belongs;
    }

    private static float dist2(float x1, float y1, float x2, float y2)
    {
        return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);
    }

    private static void extractClustersForRefinedStreamlineLocations(
            ArrayList<ArrayList<float[]>> _biggest_clusters,
            float[][][] _xy2_at_loc)
    {

        // extract biggest cluster for each branch with refined points if there is branch
        if (_xy2_at_loc!=null) {
            for (int i = 0; i < _xy2_at_loc.length; i++) { 		// loop branches
                if (_xy2_at_loc[i]!=null) { 					// there are points to cluster
                    float[] xref = _xy2_at_loc[i][0];
                    float[] yref = _xy2_at_loc[i][1];
                    _biggest_clusters.get(i).clear();
                    _biggest_clusters.set(i, largest_cluster(xref, yref, clusteringDisk)); // sets cluster as list of connected xy-disk locations
                }
                else {
                    _biggest_clusters.get(i).clear();
                }
            }
        }
        else { // all clusters are zero size
            for (int i = 0; i < _biggest_clusters.size(); i++) _biggest_clusters.get(i).clear();
        }

    }

    private static void extractLocalClusterMask(
            byte[][][] _mask,
            ArrayList<ArrayList<float[]>> _biggest_clusters,
            int _location)
    {

        // store clusters in mask at this location : local image patch 4 x 2N+1 x 2N+1
        for (int cmask = 0; cmask < _mask.length; cmask++) {
            for (int xmask = 0; xmask < _mask[0].length; xmask++) {
                for (int ymask = 0; ymask < _mask[0][0].length; ymask++) {

                    _mask[cmask][xmask][ymask] = (byte)0;

                    int xc = i2xy[_location][0];
                    int yc = i2xy[_location][1];

                    float xc_real = xc - Nimg + xmask + 0.5f; // inverse way would be: floor(real_x-Nimg)
                    float yc_real = yc - Nimg + ymask + 0.5f;
                    float rc_real = (float) 0.5;// 4 connectivity //(0.5*Math.sqrt(2));

                    // compare with extracted collection of biggest clusters
                    if (_biggest_clusters.get(cmask).size()>0) {

                        for (int j = 0; j < _biggest_clusters.get(cmask).size(); j++) {

                            float cluster_element_x = _biggest_clusters.get(cmask).get(j)[0];
                            float cluster_element_y = _biggest_clusters.get(cmask).get(j)[1];
                            float cluster_element_r = clusteringDisk;

                            boolean overlap = dist2(xc_real, yc_real, cluster_element_x, cluster_element_y) <= Math.pow(rc_real + cluster_element_r, 2);

                            if (overlap) {
                                _mask[cmask][xmask][ymask] = (byte)64;
                            }

                        }

                    }

                }

            }

        }

    }

    private static void extractRefinedStreamlineSelection(
            boolean[][] _xy2sel_at_loc,                       // output
            ArrayList<ArrayList<float[]>> _selected_xy_local, // output

            float[][][] _xy2_at_loc,
            int _loc,
            float[][][] _xy_local,

            ArrayList<ArrayList<float[]>> _biggest_clusters,
            byte[][][] _mask)
    {

        if (_xy2_at_loc!=null) {

            int root_x = i2xy[_loc][0] - Nimg;
            int root_y = i2xy[_loc][1] - Nimg;

            for (int i = 0; i < _xy2_at_loc.length; i++) {

                if (_xy2_at_loc[i]!=null) {

                    int nr_refs = _xy2_at_loc[i][0].length; // nr x coordinates

                    _xy2sel_at_loc[i] = new boolean[nr_refs];
                    Arrays.fill(_xy2sel_at_loc[i], true);

                    ArrayList<float[]> local_xy_to_add = new ArrayList<float[]>(L);

                    // which ones to throw away?
                    for (int j = 0; j < _xy2_at_loc[i][0].length; j++) { // j loops coordinates of one streamline

                        float picked_x = _xy2_at_loc[i][0][j];
                        float picked_y = _xy2_at_loc[i][1][j];
                        float picked_x_local = _xy_local[i][0][j];
                        float picked_y_local = _xy_local[i][1][j];

                        if (j>=L) {
                            _xy2sel_at_loc[i][j] = false;
                        }
                        else {

                            _xy2sel_at_loc[i][j] = true; // we expect it to be added unless...

                            // see if it belongs to the rest of the branches
                            for (int k = 0; k < _biggest_clusters.size(); k++) {
                                if (k!=i && _biggest_clusters.get(k).size()>=L/2f) {  // i - branch that we loop, k - branch that we are checking

                                    int xx = ((int) Math.floor(picked_x - root_x));
                                    xx = xx<0? 0 : xx;
                                    xx = xx>=_mask[0].length? _mask[0].length : xx;

                                    int yy = ((int) Math.floor(picked_y - root_y));
                                    yy = yy<0? 0 : yy;
                                    yy = yy>=_mask[0][0].length? _mask[0][0].length : yy;

                                    boolean belongs = (_mask[k][xx][yy]&0xff) > 0;

                                    if (belongs) _xy2sel_at_loc[i][j] = false;
                                }
                            }

                            if (_xy2sel_at_loc[i][j])
                                local_xy_to_add.add(new float[]{picked_x_local, picked_y_local});

                        }

                    }

                    _selected_xy_local.get(i).clear();
                    _selected_xy_local.set(i, local_xy_to_add);

                }
                else {
                    _xy2sel_at_loc[i] = null;
                    _selected_xy_local.get(i).clear();
                }
            }
        }
        else {
            for (int i = 0; i < _xy2sel_at_loc.length; i++) _xy2sel_at_loc[i] = null;
            for (int i = 0; i < _selected_xy_local.size(); i++) _selected_xy_local.get(i).clear();
        }

    }

    public static float getSmoothnessPercentile(int perc)
    {

        // check the distribution of (non-normalized) smoothness values
        // loop through all to count all the valid smoothness values
        int count = 0;
        float min_smoothness = Float.POSITIVE_INFINITY;
        float max_smoothness = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < smoothness.length; i++) {

            if (smoothness[i]!=null) {

                for (int j = 0; j < smoothness[i].length; j++) {

                    if (!Float.isNaN(smoothness[i][j])) {

                        count++;
                        if (smoothness[i][j]<min_smoothness) min_smoothness = smoothness[i][j];
                        if (smoothness[i][j]>max_smoothness) max_smoothness = smoothness[i][j];

                    }

                }

            }

        }

        System.out.println(count + " smoothness values found");
        System.out.println(min_smoothness + " -> min. smoothness");
        System.out.println(max_smoothness + " -> max. smoothness");

        float[] all_sth = new float[count];

        count = 0;

        for (int i = 0; i < smoothness.length; i++) {

            if (smoothness[i]!=null) {

                for (int j = 0; j < smoothness[i].length; j++) {

                    if (!Float.isNaN(smoothness[i][j])) {

                        all_sth[count++] = smoothness[i][j];

                    }

                }

            }

        }

        int t20 = get20(perc);

        return Stat.quantile(all_sth, t20, 20);

    }

    private static int get20(int percentile)
    {
        percentile = (percentile< 1)?  1 : percentile;
        percentile = (percentile>99)? 99 : percentile;
        int t20 = Math.round((20f*percentile)/100f);
        t20 = (t20< 1)? 1  : t20;
        t20 = (t20>19)? 19 : t20;
        return t20;
    }

}