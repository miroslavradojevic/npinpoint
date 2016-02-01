package com.braincadet.npinpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.PointRoi;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.util.ArrayList;

/**
 * Computation unit for 2D processing, contains precomputed values for profile extraction and peak detection
 * indexing, neighbour list (for peaks), offset list (for profiles)
 */
public class Sphere2D {

    private static float    arcRes 	        = 1.0f;
    public static float 	arcNbhood       = 2*arcRes;
    private static float    samplingStep    = 0.9f;

    public static float     T_HALF 	        = 0.50f;

    private float 	radius;
    private float   neuronDiameter;
    private float 	sigma_ratio;
    private int     N;

    private int 	limR, limT;

    private static float TWO_PI = (float) (2 * Math.PI);
    private static float ONE_PI = (float) (1 * Math.PI);

    // all discretized angle (theta) values will be indexed in a list
    // masks will contain indexes of local neighbours (necessary for peak extraction and clustering)
    // offsets used in oriented filtering are precomputed for each indexed angle
    // there is also a symmetric table that stores precomputed distance differences between each indexed direction (used in clustering methods after converging the indexes)

    public static ArrayList<Float>          theta = new ArrayList<Float>(); 	        // list of elements (theta) covering the circle
    public static ArrayList<int[]> 		    masks = new ArrayList<int[]>(); 	    // list of list indexes of the neighbours for each list element
    private static ArrayList<float[][]> 	offstXY = new ArrayList<float[][]>(); 	// list of filter offsets for each direction

    // variables for bayesian tracking
    public static ArrayList<float[]>        locsXY  = new ArrayList<float[]>();      // list of follow up location offsets (for prediction transistion)
    public static ArrayList<float[]>        vxy = new ArrayList<float[]>();          // unit vector directions corresponding to thetas

    private static float[] weights;
    public static float[][] diffs;

    /*
    *********************************************************************
     */

    public Sphere2D(float neuronDiam, float scale, float sigma_ratio)
    {

        this.radius 	= scale * neuronDiam;
        this.neuronDiameter = neuronDiam;
        this.sigma_ratio = sigma_ratio;

        this.N 	= (int) Math.ceil( ( (2 * Math.PI * radius) / arcRes) );    // N will influence theta list size, and offstXY list size
        this.limT = (int) Math.ceil(T_HALF*neuronDiameter/samplingStep);    // transversal sampling limits
        this.limR = 2 * limT + 1; // how many to take radially with given sampling step

        theta.clear();
        for (int i=0; i<N; i++) {
            theta.add( i * ( (float)(2*Math.PI) / N ) );
        }

        masks.clear();
        for (int ii=0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            // loop the rest of the elements to see if they are in the angNeighbourhood range
            ArrayList<Integer> nbrs = new ArrayList<Integer>();
            for (int jj = 0; jj<theta.size(); jj++) {

                if (jj!=ii) {

                    // check
                    float theta_test = theta.get(jj);

                    float arc_btw = arcBetweenDirections(curr_theta, theta_test);
                    if (arc_btw<=arcNbhood) {
                        nbrs.add(jj);
                    }

                }

            }

            // convert list to regular array and add
            int[] nbrsArray = new int[nbrs.size()];
            for (int iii=0; iii<nbrs.size(); iii++) {
                nbrsArray[iii] = nbrs.get(iii);
            }

            masks.add(nbrsArray);

        }

        offstXY.clear();
        weights = new float[(2*limT+1)*limR]; // (limR+1)

        float sumWgt = 0;
//		float sumWgtPos = 0;
//		float sumWgtNeg = 0;

        for (int ii = 0; ii<theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            /*
				form sampling (offset) points
			 */

            float[][] offsetsPerDirection = new float[(2*limT+1)*limR][2]; // (limR+1)

            int cnt = 0;

            for (int k=-(limR-1); k<=0; k++) {

                for (int i=-limT; i<=limT; i++) {

//                    for (int j = -limT; j<=limT; j++) {

                    float px = i * samplingStep;
//                        float py = j * samplingStep;
                    float py = k * samplingStep;

                    offsetsPerDirection[cnt][0] = px;
                    offsetsPerDirection[cnt][1] = py;
//                        offsetsPerDirection[cnt][2] = pz;

                    //*** DEFINES THE FILTER PROFILE WEIGHTS (only in first iteration (ii=0), the rest are the same)
                    if (ii==0) {
                        float dstAxis = point2line(0, 0,        0, 1,       px, py);
                        float sig = sigma_ratio * neuronDiameter;
                        weights[cnt] = (float) Math.exp(-(dstAxis*dstAxis)/(2*sig*sig));
                        sumWgt += weights[cnt];
//							weights[cnt] = (float) ((1-(dstAxis*dstAxis)/(sigm*sigm))*Math.exp(-(dstAxis*dstAxis)/(2*sigm*sigm)));
//							if (weights[cnt]>0) sumWgtPos+=weights[cnt];
//							else sumWgtNeg-=weights[cnt];
                    }

                    cnt++;

//                    }

                }

            }

            /*
				transformation for offsets before adding
			 */
            transY(radius, offsetsPerDirection);
            //rotY(-phi+HalfPI, offsetsPerDirection);
            rotZ(curr_theta, offsetsPerDirection);
            offstXY.add(offsetsPerDirection); //store

        }

        /*
            prediction locations locXY, and the directions vxy (for prior)
         */
        locsXY.clear();
        vxy.clear();
        for (int ii = 0; ii < theta.size(); ii++) {

            float curr_theta = theta.get(ii);

            int i = 0;
            int k = - (int) Math.round(limR * (1/2f));// here it scales the step

            float px = i * samplingStep;
            float py = k * samplingStep;

            float[][] locPerDirection = new float[1][2];

            locPerDirection[0][0] = px;
            locPerDirection[0][1] = py;

            transY(radius, locPerDirection);
            rotZ(curr_theta, locPerDirection);
            locsXY.add(locPerDirection[0]);

            float vx = 0;
            float vy = 0;

            float[][] vecPerDirection = new float[1][2];

            vecPerDirection[0][0] = vx;
            vecPerDirection[0][1] = vy;

            transY(1, vecPerDirection);
            rotZ(curr_theta, vecPerDirection);
            // normalize vxy (just in case)
            double vnorm = Math.sqrt( vecPerDirection[0][0]*vecPerDirection[0][0] + vecPerDirection[0][1]*vecPerDirection[0][1] );
            vecPerDirection[0][0] /= vnorm;
            vecPerDirection[0][1] /= vnorm;

            vxy.add(vecPerDirection[0]);

        }


        /*
				normalize weights
	    */
        for (int iii=0; iii<weights.length; iii++) {
            weights[iii] /= sumWgt;
//			if (weights[iii]>0)weights[iii] /= sumWgtPos;
//			else {weights[iii] /= sumWgtNeg; System.out.println("YES");}
        }

        /*
                form table with differences (used for clustering)
         */

        diffs = new float[theta.size()][theta.size()];
        for (int i = 0; i < theta.size(); i++) {
            for (int j = i; j < theta.size(); j++) {
                if (i==j) {
                    diffs[i][j] = 0;
                }
                else {
                    float theta1 = theta.get(i);
                    float theta2 = theta.get(j);
                    float dtheta = wrap_diff(theta1, theta2);
                    diffs[i][j] = radius * dtheta;
                    diffs[j][i] = radius * dtheta;
                }
            }
        }

    }

    public void getPriors(float[] _vxy, float _deg_std, float[] _store_priors)
    {

        // normalize the input direction
        double _vnorm = Math.sqrt( _vxy[0]*_vxy[0] + _vxy[1]*_vxy[1] );
        _vxy[0] /= _vnorm;
        _vxy[1] /= _vnorm;

        // take the reference direction and assign the probability to each of the transition points
        for (int i = 0; i < vxy.size(); i++) {

            float dot_prod = _vxy[0]*vxy.get(i)[0] + _vxy[1]*vxy.get(i)[1];
            dot_prod = (dot_prod>1)? 1 : dot_prod;
            double ang_rad = Math.acos(dot_prod);
            double ang_deg = ang_rad*(180f/Math.PI);

            // form gaussian weighted prior
            _store_priors[i] = (float) ((float) (1 / (_deg_std*Math.sqrt(2*Math.PI)) ) * Math.exp( -((ang_deg-0)*(ang_deg-0)) / (2*_deg_std*_deg_std) ));

        }

    }

    public ImagePlus showSampling()
    {

        int DIM = 2 * (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF,2)+Math.pow(radius, 2))) + 1;
        int CX = DIM/2;
        int CY = CX;

        ImageStack isOut = new ImageStack(DIM, DIM);
        Overlay ov = new Overlay();

        for (int i=0; i<offstXY.size(); i++) {
            isOut.addSlice(new ByteProcessor(DIM, DIM));

            // center
            OvalRoi p = new OvalRoi(CX+0.5 - .5, CY+0.5 -.5, 1, 1);
            p.setPosition(i+1);
            p.setFillColor(Color.RED);
            p.setFillColor(Color.RED);
            ov.add(p);

            // sampling
            for (int i1=0; i1<offstXY.get(i).length; i1++) {

                float offX = offstXY.get(i)[i1][0];
                float offY = offstXY.get(i)[i1][1];

                PointRoi p1 = new PointRoi(CX+offX+.5, CY+offY+.5);
                p1.setPosition(i+1);
                ov.add(p1);

            }

        }

        ImagePlus outIm = new ImagePlus("offsets", isOut);
        outIm.setOverlay(ov);
        return outIm;

    }

    public void exportSampling(String file_path)
    {

        FileSaver fs = new FileSaver(showSampling());
        fs.saveAsTiffStack(file_path);

    }

    public ImagePlus showWeights()
    {

        float sum = 0;
        for (int k=0; k<weights.length; k++) sum += weights[k];
        float sig = sigma_ratio * neuronDiameter;
        return new ImagePlus("weights,sigma="+ IJ.d2s(sig, 2)+",sum="+IJ.d2s(sum,2), new FloatProcessor((2*limT+1), limR, weights)); // (limR+1)

    }

    public void exportWeights(String file_path)
    {

        FileSaver fs = new FileSaver(showWeights());
        fs.saveAsTiff(file_path);

    }

    public int getProfileLength()
    {
        return offstXY.size();
    }

    public int getOuterRadius()
    {
        return (int) Math.ceil(Math.sqrt(Math.pow(neuronDiameter*T_HALF, 2)+Math.pow(radius, 2)));
    }

    public float getScale()
    {
        return radius / neuronDiameter;
    }

    public float getInitialRadialSeparation()
    {
        return radius - neuronDiameter;
    }

    public float getNeuronDiameter()
    {
        return neuronDiameter;
    }

    public float getSphereRadius()
    {
        return radius;
    }

    public short extractProfile(int profileIdx, float atX, float atY, float[][] _inimg_xy)
    {
        // one element filter output (indexed with profileIdx)

        float value = 0;

        for (int offsetIdx=0; offsetIdx<offstXY.get(profileIdx).length; offsetIdx++) {

            float x_offs_pix = offstXY.get(profileIdx)[offsetIdx][0];
            float y_offs_pix = offstXY.get(profileIdx)[offsetIdx][1];

            float imgVal = Interpolator.interpolateAt(atX + x_offs_pix, atY + y_offs_pix, _inimg_xy); // , atZ + z_offs_lay
            value += weights[offsetIdx] * imgVal;

        }

        return  (short) ((int) ((value/255f)*65535f)); // &  0xffff

    }

//	public void peakCoords_4xXY(short[] _profile,                           	// main input
//                                int[] start_pts, int[] end_pts,          		// auxiliary. arrays (to avoid allocating them inside method each time) - end_pts is also an output
//                                int atX, int atY,                           	// for global coordinate outputs
//                                float[][] _inimg_xy, int[][] _xy2i,         	// to calculate median() along line
//                                int[][] peaks_loc_xy, float[][] peaks_ang_theta)
//	{	// outputs
//
//        // this is the block that extracts 4 strongest peaks in global 2d coordinates (that's why atX, atY at the input) and profile indexes
//        // ranking can be based on median between locations of number of points that converged to each cluster
//        // in order to rank using medians - image array, look-up tables are added to be able to compute medians
//        // otherwise ranking can be based on number of convergence points - then all the image data is not necessary (position is enough)
//        // this implementation will have the number of convergence points but use median along the line as ranking estimate finally
//        runLineSearch(
//                start_pts,
//                _profile,
//                MAX_ITER,
//                EPSILON,
//                end_pts);
//
//        // cluster end_pts together
//		int[] labs = clustering(end_pts, diffs, arcNbhood);
//
//		// extract the cluster centroids out  -> <theta, weight>
//		ArrayList<float[]> clss = extracting(labs, end_pts, theta);
//
//		appendClusters(atX, atY, _xy2i, _inimg_xy, clss, peaks_loc_xy, peaks_ang_theta);
//
//	}

    /*
    *********************************************************************
     */

//	private void appendClusters(
//									   int atX, int atY,
//									   int[][] 		_xy2i,  						// lookup table
//									   float[][] 	_inimg_xy,                     	// input image
//									   ArrayList<float[]> clusters_to_append,       // <theta, weight>
//									   int[][] 		destination_locs,               // 4x2
//									   float[][] 	destination_angs)               // 4x1
//	{
//
//        // _inimg_xy was used just to be input for median along the line anf to
//        // be aware of the image dimensions - check if the appended location is within the image
//        int W = _inimg_xy.length;
//        int H = _inimg_xy[0].length;
//
//		// take top 4 and store them according to the importance (criteria: median along the connecting line, or convergence iterations)
//		float[] medAlongLin = new float[destination_locs.length];
//		Arrays.fill(medAlongLin, -1f);
//
//		for (int t=0; t<clusters_to_append.size(); t++) { // check every peak theta angle
//
//			float peak_theta = clusters_to_append.get(t)[0];   // value in [rad] theta !!!!!!!!!
//            int[] currentPeaksXY   = new int[2];
//
//			float x_peak_pix = atX + getX(peak_theta);
//			float y_peak_pix = atY + getY(peak_theta);
//
//            boolean modified = false;
//            float maxWeight = Float.MIN_VALUE;
//
//            if (CHECK_NBRS_WHEN_ASSIGNING_PEAK_LOCS) {
//
//                /*
//				pick the best out of 4 neighbours (compensate location inaccuracy), can be nothing if it's either in background or out of the image
//			    */
//
//                int x_peak_pix_base = (int) Math.floor(x_peak_pix);
//                int y_peak_pix_base = (int) Math.floor(y_peak_pix);
//
//                for (int ii = 0; ii <= 1; ii++) { // check around peak
//                    for (int jj = 0; jj <= 1; jj++) {
//
//                        int check_x = x_peak_pix_base + ii;
//                        int check_y = y_peak_pix_base + jj;
//
//                        if (check_x>=0 && check_x<W && check_y>=0 && check_y<H) { // is in image
//
//                            float currMedian = PeakAnalyzer2D.medianAlongLine(	atX, atY, check_x, check_y, _inimg_xy );
//
//                            if (_xy2i[check_x][check_y]!=-1)  {  // ensures retrieval from the list
//
//                                if (currMedian>maxWeight) {
//                                    // set this one as peak
//                                    currentPeaksXY[0] = check_x;
//                                    currentPeaksXY[1] = check_y;
//
//                                    modified = true;
//
//                                    // update maxMedian
//                                    maxWeight = currMedian;
//                                }
//
//                            }
//                        }
//                    }
//                }
//
//            }
//            else {
//
//                // rounded locations, don't check the neighbours
//                int x_peak_pix_rounded = (int) Math.round(x_peak_pix);
//                int y_peak_pix_rounded = (int) Math.round(y_peak_pix);
//
//                if (x_peak_pix_rounded>=0 && x_peak_pix_rounded<W && y_peak_pix_rounded>=0 && y_peak_pix_rounded<H) {
//                    float currMedian = PeakAnalyzer2D.medianAlongLine(	atX, atY, x_peak_pix_rounded, y_peak_pix_rounded, _inimg_xy );
//                    if (_xy2i[x_peak_pix_rounded][y_peak_pix_rounded]!=-1) {
//                        if (currMedian>maxWeight) { // dummy thing, just to be consistent
//                            currentPeaksXY[0] = x_peak_pix_rounded;
//                            currentPeaksXY[1] = y_peak_pix_rounded;
//                            modified = true;
//                            maxWeight = currMedian;
//                        }
//                    }
//                }
//
//            }
//
//			if (modified) {
//
//			    // if there was a result append it to sorted list of peaks around this location
//				// insert currentPeaksXY and peak_theta -> to the list of 4 (destination_locs, destination_angs)
//				for (int k = 0; k < 4; k++) {
//
//					if (medAlongLin[k] == -1f) {            			// store immediately
//
//						destination_locs[k][0]    = currentPeaksXY[0];
//						destination_locs[k][1]    = currentPeaksXY[1];
//
//						medAlongLin[k]  = maxWeight;
//
//						destination_angs[k][0]     = peak_theta;
//
//						break;
//
//					}
//					else if (maxWeight > medAlongLin[k]) {
//
//						// shift the rest first
//						for (int kk = 4-2; kk>=k; kk--) {   // shift them from the one before last, shift back (last one dissapears)
//
//							destination_locs[kk+1][0] = destination_locs[kk][0];
//							destination_locs[kk+1][1] = destination_locs[kk][1];
//
//							medAlongLin[kk+1] = medAlongLin[kk];
//
//							destination_angs[kk+1][0] = destination_angs[kk][0];
//
//						}
//
//						// store it at k
//						destination_locs[k][0] = currentPeaksXY[0];
//						destination_locs[k][1] = currentPeaksXY[1];
//
//						medAlongLin[k] = maxWeight;
//
//						destination_angs[k][0] = peak_theta;
//
//						break;
//
//					}
//					else {
//
//						// if smaller, loop further
//
//					}
//
//				}
//
//			} // modified
//
//		}
//
//
//	}

    private void rotZ(float ang, float[][] coords) {
        for (int i=0; i<coords.length; i++) {
            float x_temp = coords[i][0];
            float y_temp = coords[i][1];
            coords[i][0] = x_temp * (float) Math.cos(ang) - y_temp * (float) Math.sin(ang);
            coords[i][1] = x_temp * (float) Math.sin(ang) + y_temp * (float) Math.cos(ang);
        }
    }

    private void transY(float ty, float[][] coords){
        for (int i=0; i<coords.length; i++){
            coords[i][1] += ty;
        }
    }

    private float point2line(float n1x, float n1y,  // float n1z,
                             float n2x, float n2y,  // float n2z,
                             float px,  float py    //, float pz
    )
    {

        float d = 0;

        double[] p_b = new double[2];

        //double[] n21 = new double[3];
        float n21Len = (float) Math.sqrt(Math.pow(n2x-n1x,2)+Math.pow(n2y-n1y,2)); // +Math.pow(n2z-n1z,2)
        float n21x = (n2x-n1x)/n21Len;
        float n21y = (n2y-n1y)/n21Len;
//        float n21z = (n2z-n1z)/n21Len;

        float proj = (px - n1x) * n21x + (py - n1y) * n21y; // + (pz - n1z) * n21z; // dot prod

        p_b[0] = -(px - n1x) + proj * n21x;
        p_b[1] = -(py - n1y) + proj * n21y;
//        p_b[2] = -(pz - n1z) + proj * n21z;

        return (float) Math.sqrt(p_b[0]*p_b[0] + p_b[1]*p_b[1]); // + p_b[2]*p_b[2]

    }

    private float arcBetweenDirections(float theta1, float theta2)
    {

        float x1 = getX(theta1);
        float y1 = getY(theta1);

        float x2 = getX(theta2);
        float y2 = getY(theta2);

        return radius * (float) Math.acos( (x1*x2+y1*y2)/(radius * radius) );

    }

//    private static float getX(float r, float theta){return (-1) * r * (float) Math.sin(theta);}
//
//    private static float getY(float r, float theta){return (+1) * r * (float) Math.cos(theta);}


    // getX() and getY() are used to correlate angle with spatial coordinates
    public float getX(float theta) {return (-1) * radius * (float) Math.sin(theta);}

    public float getY(float theta) {return (+1) * radius * (float) Math.cos(theta);}

    public static ArrayList<float[]> extracting(int[] labels, int[] idxs, ArrayList<Float> vals) {

        // loop obtained labels (labels & idxs should have the same length)

        boolean[] checked = new boolean[idxs.length];      // auxiliary
        ArrayList<float[]> out = new ArrayList<float[]>(); // allocate the output

        for (int i = 0; i < idxs.length; i++) {
            if (!checked[i]) {
                // this is the first value
                float centroid  = vals.get(idxs[i]); // vals[ idxs[i] ];
                float shifts    = 0; // naturally 0 shift for the first one
                int count = 1;
                checked[i] = true;

                // check the rest
                for (int j = i+1; j < idxs.length; j++) {
                    if (!checked[j]) {
                        if (labels[j]==labels[i]) {

                            // clustering said they're together
                            // important that vals and centroid values are all wrapped in [0, 2PI) range
                            float add_value = vals.get(idxs[j]);
                            float add_diff = wrap_diff(centroid, add_value);
                            if (centroid<ONE_PI) {
                                // there is always pi values on right
                                if (add_value>=centroid & add_value<centroid+ONE_PI) {
                                    // sign is (+)
                                }
                                else {
                                    // sign is (-)
                                    add_diff = (-1) * add_diff;
                                }

                            }
                            else {
                                // >= ONE_PI
                                // there is always pi values on left
                                if (add_value<=centroid & add_value>centroid-ONE_PI) {
                                    // sign is (-)
                                    add_diff = (-1) * add_diff;
                                }
                                else {
                                    // sign in (+)
                                }
                            }

                            shifts += add_diff;
                            count++;
                            checked[j] = true;

                        }
                    }
                }

                centroid += shifts/count;
//                centroid = ;
                out.add(new float[]{wrap_0_2PI(centroid), count});  // outputs centroid in [rad]

            }
        }

        return out;

    }

    private static float wrap_diff(float theta_1, float theta_2) { // wraps the angle difference theta_1, theta_2 range [0, 2PI)

        float d = theta_1 - theta_2;

        d = (d>0)? d : (-d) ;

        d = (d>Math.PI)? (float) (2*Math.PI-d) : d ;

        return d;

    }

    private static float wrap_0_2PI(float ang) {
        float out = ang;
        while (out<0) {
            out += TWO_PI;
        }
        while (out>=TWO_PI) {
            out -= TWO_PI;
        }
        return out;
    }

//    private static int runOneMax(int curr_pos, short[] _input_profile) {
//        int 	    new_pos     = curr_pos;
//        int 		max	 		= Integer.MIN_VALUE;
//
//        // curr_pos will define the set of neighbouring indexes
//        int[] neighbour_idxs = masks.get(curr_pos);
//
//        for (int i=0; i<neighbour_idxs.length; i++) {
//            int neighbour_idx = neighbour_idxs[i];
//            int read_value = (int) (_input_profile[neighbour_idx] & 0xffff);
//            if (read_value>max) {
//                max = read_value;
//                new_pos = neighbour_idx;
//            }
//        }
//
//        return new_pos;
//    }

//    private static void runLineSearch(
//            int[] 	    start_idxs,
//            short[] 	input_profile,
//            int 		max_iter,
//            int 	    epsilon,
//            int[] 	    end_idxs // same length as start (this would be output)
//    )
//    {
//
//        // initialize output array
//        for (int i = 0; i < end_idxs.length; i++) {
//            end_idxs[i] = start_idxs[i];
//        }
//
//        // iterate each of the elements of the output
//        for (int i = 0; i < end_idxs.length; i++) {
//
//            int iter = 0;
//            int d;
//
//            do{
//
//                int new_pos = runOneMax(end_idxs[i], input_profile);
//                int pre_value = input_profile[end_idxs[i]] & 0xffff;
//                int new_value = input_profile[new_pos]     & 0xffff;
//                d = Math.abs(new_value - pre_value);
//                end_idxs[i] = new_pos;
//                iter++;
//            }
//            while(iter < max_iter && d > epsilon);
//
//        }
//
//    }

//    private static int[] clustering(int[] idxs, float[][] dists, float threshold_dists)   // essentially clustering of the indexes
//    {
//        // indxs represent indexes of values that need to be clustered
//        // intended to place here indexes after the convergence
//        // dists are the distances
//        // threshold_dists is the distance limit
//        // output is list of unique labels
//
//        int[] labels = new int[idxs.length];
//        for (int i = 0; i < labels.length; i++) labels[i] = i;  // initialize the output
//
//        //System.out.println("INIT. LABELS:");
//        //System.out.println(Arrays.toString(labels));
//
//        for (int i = 0; i < idxs.length; i++) {
//
//            // one versus the rest
//            for (int j = 0; j < idxs.length; j++) {
//
//                // check the rest of the values
//                if (i != j) {
//
//                    int idx_i = idxs[i]; // will be used to read diff from the table
//                    int idx_j = idxs[j]; //
//
//                    if (dists[idx_i][idx_j]<=threshold_dists) {
//
//                        if (labels[j] != labels[i]) {
//                            // propagate the label
//                            int currLabel = labels[j];
//                            int newLabel  = labels[i];
//
//                            labels[j] = newLabel;
//
//                            //set all that also were currLabel to newLabel
//                            for (int k = 0; k < labels.length; k++)
//                                if (labels[k]==currLabel)
//                                    labels[k] = newLabel;
//
//                        }
//
//                    }
//
//                }
//
//            }
//
//        }
//
//        //System.out.println("OUT LABELS:");
////		for (int ii = 0; ii < labels.length; ii++)
////			System.out.print(labels[ii]+" ");
//        //System.out.println(Arrays.toString(labels));
//
//        return labels;
//
//    }

}
