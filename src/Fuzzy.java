import auxiliary.Stat;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.*;
import java.util.Arrays;

public class Fuzzy {

    // 3 inputs:
    // lhood		(up to 4 peaks of the filtering profile: template convolved to obtain profile values)
    // smoothness   (integral of the second derivative of the refined locations)
    // ncc 			(ncc of the template used with the underlying image contents)
    // 2 outputs:
    // output 1 - saying whether other streamline (that takes lhood, smoothness and ncc) is ON, NONE or OFF
    // output 2 - saying whether location is END, NONE, or BIF (up to three streamlines taken)

    // ncc range (values below 0 are cut off)
    float ncc_start = 0;
    float ncc_end 	= 1;
    float ncc_high; // trapezoid
    float ncc_low;

    // likelihood range (for visualizaiton)
    float likelihood_start 	= 0;
    float likelihood_end 	= 1;
    float likelihood_high; // trapezoid
    float likelihood_low;

    // smoothness range (for visualization)
    float smoothness_start = 0;
    float smoothness_end = Float.NaN; // because it will be inferred from the input arguments, needs to be initiated
    float smoothness_high; // will be given as an argument
    float smoothness_low;  // argument

    float 	output_sigma;

    // output 1 can be off, none or on
    float 	output1_range_start;
    float 	output1_range_ended;
    float 	output1_off  = 0f;
    float 	output1_none = 1f;
    float   output1_on 	 = 2f;

    // output 3 can be end, none or jun
    float 	output2_range_start;
    float 	output2_range_ended;
    float 	output2_endpoint = 1f;
    float   output2_nonpoint = 2f;
    float   output2_bifpoint = 3f;

    int      N1 = 100;         		// number of points to aggregate for each output, hardcoded, need to cover the range (x_min, x_max) well
    int		 N2 = N1;

    // fuzzification parameters
    float[] 	x1;      	    	// serve as x axis for membership, aggregation functions for output1 (range output1_range_start to output1_range_ended)
    float[] 	x2;      	    	// serve as x axis for membership, aggregation functions for output2 (range output2_range_start to output2_range_ended)
    float x_min = Float.POSITIVE_INFINITY, x_max = Float.NEGATIVE_INFINITY;

    public float[] 	agg1;            // serve as aggregation for output1 (off, none, on)
    public float[] 	agg2;            // serve as aggregation for output2 (end, none, bif)

    // output categories
    float[] 	v_off;
    float[]		v_none;
    float[] 	v_on;

    float[]		v_endpoint;
    float[]		v_nonpoint;
    float[]		v_bifpoint;

    // membership fuzzy sets for output
    float[] 	q_off;
    float[] 	q_none;
    float[] 	q_on;
    float[] 	q_endpoint;
    float[] 	q_nonpoint;
    float[] 	q_bifpoint;

    // detection log
    public boolean verbose = false; // if verbose, the stack will be appended with score for each detection afterwards
    private static int plotw = new Plot("","","",new float[1],new float[1]).getSize().width;
    private static int ploth = new Plot("","","",new float[1],new float[1]).getSize().height;
    public ImageStack fls_steps = new ImageStack(plotw, ploth); // will append current decisions each time they're done if we decide to set verbose=true

    public Fuzzy(
            float _ncc_high,
            float _ncc_low,
            float _likelihood_high,
            float _likelihood_low,
            float _smoothness_high,
            float _smoothness_low,
            float _output_sigma
    )
    {
        ncc_high 	= _ncc_high;
        ncc_low 	= _ncc_low;
        likelihood_high = _likelihood_high;
        likelihood_low 	= _likelihood_low;
        smoothness_high = _smoothness_high; // low values mean high smoothness
        smoothness_low = _smoothness_low;
        output_sigma = _output_sigma;

        smoothness_end = smoothness_low + (smoothness_low-smoothness_high)/2;

        output1_range_start = output1_off-3*output_sigma;
        output1_range_ended = output1_on+3*output_sigma;

        output2_range_start = output2_endpoint-3*output_sigma;
        output2_range_ended = output2_bifpoint+3*output_sigma;

        x1 = new float[N1]; // output1 range
        x2 = new float[N2]; // output2 range

        for (int i=0; i<N1; i++) x1[i] = output1_range_start + i * ( (output1_range_ended-output1_range_start) / (N1-1));
        for (int i=0; i<N2; i++) x2[i] = output2_range_start + i * ( (output2_range_ended-output2_range_start) / (N2-1));

        for (int i=0; i<N1; i++) {if(x1[i]<x_min) x_min = x1[i]; if(x1[i]>x_max) x_max = x1[i];}
        for (int i=0; i<N2; i++) {if(x2[i]<x_min) x_min = x2[i]; if(x2[i]>x_max) x_max = x2[i];}

        // output membership functions
        q_off = new float[N1];
        for (int i=0; i<N1; i++)
            q_off[i] = (float) Math.exp(-Math.pow(x1[i] - output1_off, 2) / (2 * Math.pow(output_sigma, 2)));

        q_none = new float[N1];
        for (int i=0; i<N1; i++)
            q_none[i] = (float) Math.exp(-Math.pow(x1[i] - output1_none, 2) / (2 * Math.pow(output_sigma, 2)));

        q_on = new float[N1];
        for (int i=0; i<N1; i++)
            q_on[i] = (float) Math.exp(-Math.pow(x1[i] - output1_on, 2) / (2 * Math.pow(output_sigma, 2)));

        q_endpoint = new float[N2];
        for (int i=0; i<N2; i++)
            q_endpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));

        q_nonpoint = new float[N2];
        for (int i=0; i<N2; i++)
            q_nonpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));

        q_bifpoint = new float[N2];
        for (int i=0; i<N2; i++)
            q_bifpoint[i] = (float) Math.exp(-Math.pow(x2[i] - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));

        // auxiliary arrays used for calculations
        agg1 = new float[N1];
        agg2 = new float[N2];

        v_off = new float[N1];
        v_none = new float[N1];
        v_on  = new float[N1];

        v_endpoint = new float[N2];
        v_nonpoint = new float[N2];
        v_bifpoint = new float[N2];

        verbose = false;
        // clear the log
        for (int i = 0; i < fls_steps.getSize(); i++) fls_steps.deleteSlice(i + 1);

    }

    public void clearLog()
    {
//		System.out.println("size before " + fls_steps.getSize());
        while (fls_steps.getSize()>0) {
            fls_steps.deleteLastSlice();
        }
//		for (int i = 0; i < fls_steps.getSize(); i++) {fls_steps.deleteSlice(i + 1);
//			System.out.println(""+(i+1));}
    }

    /*
    fuzzification of the ncc input
     */
    private float h_ncc_high(float ncc)
    {
        if (ncc>=ncc_high)
            return 1;
        else if (ncc<ncc_high && ncc>=ncc_low)
            return (ncc-ncc_low)/(ncc_high-ncc_low);
        else
            return 0;
    }

    private float h_ncc_low(float ncc)
    {
        if (ncc<ncc_low)
            return 1;
        else if (ncc>=ncc_low && ncc<ncc_high)
            return (ncc-ncc_high)/(ncc_low-ncc_high);
        else
            return 0;
    }

    /*
    fuzzification of the likelihood input
     */
    private float h_likelihood_high(float likelihood)
    {
        if (likelihood>=likelihood_high)
            return 1;
        else if (likelihood<likelihood_high && likelihood>=likelihood_low)
            return (likelihood-likelihood_low)/(likelihood_high-likelihood_low);
        else
            return 0;
    }

    private float h_likelihood_low(float likelihood)
    {
        if (likelihood<likelihood_low)
            return 1;
        else if (likelihood>=likelihood_low && likelihood<likelihood_high)
            return (likelihood-likelihood_high)/(likelihood_low-likelihood_high);
        else
            return 0;
    }

    /*
    fuzzification of the smooothness input
     */
    private float h_smoothness_high(float smoothness)
    {
        if (smoothness<=smoothness_high)
            return 1;
        else if (smoothness>smoothness_high && smoothness<=smoothness_low)
            return (smoothness_low-smoothness)/(smoothness_low-smoothness_high);
        else
            return 0;
    }

    private float h_smoothness_low(float smoothness)
    {
        if (smoothness>smoothness_low)
            return 1;
        else if (smoothness<=smoothness_low && smoothness>smoothness_high)
            return (smoothness-smoothness_high)/(smoothness_low-smoothness_high);
        else
            return 0;
    }

    private float[] fi_off(float z)
    {
        Arrays.fill(v_off, 0);

        for (int i=0; i<N1; i++)
            if (q_off[i]<=z)
                v_off[i] = q_off[i];
            else
                v_off[i] = z;

        return v_off;
    }

    private float[] fi_none(float z)
    {
        Arrays.fill(v_none, 0);

        for (int i=0; i<N2; i++)
            if (q_none[i]<=z)
                v_none[i] = q_none[i];
            else
                v_none[i] = z;

        return v_none;
    }

    private float[] fi_on(float z)
    {
        Arrays.fill(v_on, 0);

        for (int i=0; i<N1; i++)
            if (q_on[i]<=z)
                v_on[i] = q_on[i];
            else
                v_on[i] = z;

        return v_on;
    }

    private float[] fi_endpoint(float z)
    {
        Arrays.fill(v_endpoint, 0);

        for (int i=0; i<N2; i++)
            if (q_endpoint[i]<=z)
                v_endpoint[i] = q_endpoint[i];
            else
                v_endpoint[i] = z;

        return v_endpoint;
    }

    private float[] fi_nonpoint(float z)
    {
        Arrays.fill(v_nonpoint, 0);

        for (int i=0; i<N2; i++)
            if (q_nonpoint[i]<=z)
                v_nonpoint[i] = q_nonpoint[i];
            else
                v_nonpoint[i] = z;
        return v_nonpoint;
    }

    private float[] fi_bifpoint(float z)
    {
        Arrays.fill(v_bifpoint, 0);

        for (int i=0; i<N2; i++)
            if (q_bifpoint[i]<=z)
                v_bifpoint[i] = q_bifpoint[i];
            else
                v_bifpoint[i] = z;

        return v_bifpoint;
    }

    /*
    *  pointer branch
    *  corresponds to the branch with highest likelihood 1
    *  each method has two versions - one where it gives defuzzified value
    *  and one where it separates the membership to output categories
     */

    /****************************************************************************/

    private float branchStrength(float ncc_in, float likelihood_in, float smoothness_in)
    {
        Arrays.fill(agg1, 0);
        float[] cur;
        float mu;

        mu = min(h_likelihood_low(likelihood_in),   h_smoothness_low(smoothness_in),    h_ncc_low(ncc_in));	    cur = fi_off(mu);	accumulate(cur, agg1);
        mu = min(h_likelihood_low(likelihood_in),   h_smoothness_low(smoothness_in),    h_ncc_high(ncc_in));	cur = fi_off(mu);	accumulate(cur, agg1);
        mu = min(h_likelihood_low(likelihood_in),   h_smoothness_high(smoothness_in),   h_ncc_low(ncc_in));	    cur = fi_off(mu);	accumulate(cur, agg1);
//        mu = min(h_likelihood_low(likelihood_in),   h_smoothness_high(smoothness_in),   h_ncc_high(ncc_in));	cur = fi_none(mu);	accumulate(cur, agg1); //?

        mu = min(h_likelihood_high(likelihood_in),  h_smoothness_low(smoothness_in),    h_ncc_low(ncc_in));	    cur = fi_none(mu);	accumulate(cur, agg1);
        mu = min(h_likelihood_high(likelihood_in),  h_smoothness_low(smoothness_in),    h_ncc_high(ncc_in));	cur = fi_none(mu);	accumulate(cur, agg1);
        mu = min(h_likelihood_high(likelihood_in),  h_smoothness_high(smoothness_in),   h_ncc_low(ncc_in));	    cur = fi_none(mu);	accumulate(cur, agg1);

        mu = min(h_likelihood_high(likelihood_in),  h_smoothness_high(smoothness_in),   h_ncc_high(ncc_in));	cur = fi_on(mu);	accumulate(cur, agg1);

        return get_agg1_centroid();
    }

    public void branchStrength(float ncc_in, float likelihood_in, float smoothness_in, float[] output_2_OFF_NONE_ON) // output_off_none_on: [off_score, none_score, on_score]
    {



        float defuzz = branchStrength(ncc_in, likelihood_in, smoothness_in);   // will update agg1
        output_2_OFF_NONE_ON[0] = (float) Math.exp(-Math.pow(defuzz - output1_off, 2)   / (2 * Math.pow(output_sigma, 2)));
        output_2_OFF_NONE_ON[1] = (float) Math.exp(-Math.pow(defuzz - output1_none, 2)  / (2 * Math.pow(output_sigma, 2)));
        output_2_OFF_NONE_ON[2] = (float) Math.exp(-Math.pow(defuzz - output1_on, 2)    / (2 * Math.pow(output_sigma, 2)));

        String title = 	"OFF=" + IJ.d2s(output_2_OFF_NONE_ON[0], 2) +
                "?="   + IJ.d2s(output_2_OFF_NONE_ON[1], 2) +
                "ON="  + IJ.d2s(output_2_OFF_NONE_ON[2], 2);
        if (verbose)  appendAgg1(defuzz, title);
    }

    /****************************************************************************/
    public void critpointScore(float ncc_1, float likelihood_1, float smoothness_1,
                               float[] output_3_END_NON_JUN,
                               float[] branch_saliency)
    {

        if (verbose) appendFuzzification(ncc_1, likelihood_1, smoothness_1);

        float[] t = new float[3];
        branchStrength(ncc_1, likelihood_1, smoothness_1, t);   // will update agg1 & append if verbose is true
        float b1_is_off 	= t[0];
        float b1_is_none	= t[1];
        float b1_is_on 		= t[2];
        branch_saliency[0]  = b1_is_on;

        branch_saliency[1]  = Float.NaN;
        branch_saliency[2]  = Float.NaN;
        branch_saliency[3]  = Float.NaN;

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        mu = b1_is_off;     cur = fi_nonpoint(mu);   	accumulate(cur, agg2);
        mu = b1_is_none;    cur = fi_nonpoint(mu);   	accumulate(cur, agg2);
        mu = b1_is_on;      cur = fi_endpoint(mu);    	accumulate(cur, agg2);

        float defuzz = get_agg2_centroid();

        output_3_END_NON_JUN[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_3_END_NON_JUN[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_3_END_NON_JUN[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));


        // additional tweak - don't allow it to be higher than agg2 value at the defuzz point
        float limit = get_agg2(defuzz);

        output_3_END_NON_JUN[0] = Math.min(limit, output_3_END_NON_JUN[0]);
        output_3_END_NON_JUN[1] = Math.min(limit, output_3_END_NON_JUN[1]);
        output_3_END_NON_JUN[2] = Math.min(limit, output_3_END_NON_JUN[2]);

        String title =  "END=" + IJ.d2s(output_3_END_NON_JUN[0],2) +
                "NON=" + IJ.d2s(output_3_END_NON_JUN[1],2) +
                "JUN=" + IJ.d2s(output_3_END_NON_JUN[2],2);
        if (verbose) appendAgg2(defuzz, title);

    }
    /****************************************************************************/
    public void critpointScore(
            float ncc_1, float likelihood_1, float smoothness_1,
            float ncc_2, float likelihood_2, float smoothness_2,
            float[] output_endpoint_nonpoint_bifpoint,
            float[] branch_saliency) // branch_saliency has to have length 4
    {

        if (verbose) appendFuzzification(ncc_1, likelihood_1, smoothness_1);

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, smoothness_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];
        branch_saliency[0] = b1_is_on;

        if (verbose) appendFuzzification(ncc_2, likelihood_2, smoothness_2);

        branchStrength(ncc_2, likelihood_2, smoothness_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];
        branch_saliency[1] 	= b2_is_on;

        branch_saliency[2] 	= Float.NaN;
        branch_saliency[3] 	= Float.NaN;

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        mu = min(b1_is_off, 	b2_is_off);     cur = fi_nonpoint(mu); 	accumulate(cur, agg2);
        mu = min(b1_is_on, 		b2_is_on);      cur = fi_nonpoint(mu); 	accumulate(cur, agg2);

        mu = min(b1_is_off, 	b2_is_on);      cur = fi_endpoint(mu); 	accumulate(cur, agg2);
        mu = min(b1_is_on, 		b2_is_off);     cur = fi_endpoint(mu); 	accumulate(cur, agg2);

        mu = max(b1_is_none, 	b2_is_none);	cur = fi_nonpoint(mu); 	accumulate(cur, agg2);

        float defuzz =  get_agg2_centroid();

        output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));



        // additional tweak - don't allow it to be higher than agg2 value at the defuzz point
        float limit = get_agg2(defuzz);

        output_endpoint_nonpoint_bifpoint[0] = Math.min(limit, output_endpoint_nonpoint_bifpoint[0]);
        output_endpoint_nonpoint_bifpoint[1] = Math.min(limit, output_endpoint_nonpoint_bifpoint[1]);
        output_endpoint_nonpoint_bifpoint[2] = Math.min(limit, output_endpoint_nonpoint_bifpoint[2]);




        String title = 			"END=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[0],2) +
                "NON=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[1],2) +
                "JUN=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[2],2);

        if (verbose) appendAgg2(defuzz, title);

    }
    /****************************************************************************/
    public void critpointScore(
            float ncc_1, float likelihood_1, float smoothness_1,
            float ncc_2, float likelihood_2, float smoothness_2,
            float ncc_3, float likelihood_3, float smoothness_3,
            float[] output_endpoint_nonpoint_bifpoint,
            float[] branch_saliency)
    {

        if (verbose) appendFuzzification(ncc_1, likelihood_1, smoothness_1);

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, smoothness_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];
        branch_saliency[0]  = b1_is_on;

        if (verbose) appendFuzzification(ncc_2, likelihood_2, smoothness_2);

        branchStrength(ncc_2, likelihood_2, smoothness_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];
        branch_saliency[1]  = b2_is_on;

        if (verbose) appendFuzzification(ncc_3, likelihood_3, smoothness_3);

        branchStrength(ncc_3, likelihood_3, smoothness_3, t);
        float b3_is_off 	= t[0];
        float b3_is_none 	= t[1];
        float b3_is_on 		= t[2];
        branch_saliency[2]  = b3_is_on;

        branch_saliency[3] 	= Float.NaN;

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        mu = min(b1_is_off, b2_is_off, b3_is_off);     	cur = fi_nonpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_off, b2_is_on,  b3_is_on);      	cur = fi_nonpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, 	b2_is_off, b3_is_on);       cur = fi_nonpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, 	b2_is_on,  b3_is_off);      cur = fi_nonpoint(mu); accumulate(cur, agg2);

        mu = min(b1_is_off, b2_is_off, b3_is_on);      	cur = fi_endpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_off, b2_is_on,  b3_is_off);     	cur = fi_endpoint(mu); accumulate(cur, agg2);
        mu = min(b1_is_on, 	b2_is_off, b3_is_off);      cur = fi_endpoint(mu); accumulate(cur, agg2);

        mu = min(b1_is_on, 	b2_is_on,  b3_is_on);       cur = fi_bifpoint(mu); accumulate(cur, agg2);

        mu = max(b1_is_none, b2_is_none, b3_is_none);  	cur = fi_nonpoint(mu); accumulate(cur, agg2);

        float defuzz = get_agg2_centroid();

        output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));


        // additional tweak - don't allow it to be higher than agg2 value at the defuzz point
        float limit = get_agg2(defuzz);

        output_endpoint_nonpoint_bifpoint[0] = Math.min(limit, output_endpoint_nonpoint_bifpoint[0]);
        output_endpoint_nonpoint_bifpoint[1] = Math.min(limit, output_endpoint_nonpoint_bifpoint[1]);
        output_endpoint_nonpoint_bifpoint[2] = Math.min(limit, output_endpoint_nonpoint_bifpoint[2]);




        String title =  "END=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[0],2) +
                "NON=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[1],2) +
                "JUN=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[2],2);
        if (verbose) appendAgg2(defuzz, title);

    }
    /****************************************************************************/
    public void critpointScore(
            float ncc_1, float likelihood_1, float smoothness_1,
            float ncc_2, float likelihood_2, float smoothness_2,
            float ncc_3, float likelihood_3, float smoothness_3,
            float ncc_4, float likelihood_4, float smoothness_4,
            float[] output_endpoint_nonpoint_bifpoint,
            float[] branch_saliency)
    {

        if (verbose) appendFuzzification(ncc_1, likelihood_1, smoothness_1);

        float[] t = new float[3];

        branchStrength(ncc_1, likelihood_1, smoothness_1, t);
        float b1_is_off 	= t[0];
        float b1_is_none 	= t[1];
        float b1_is_on 		= t[2];
        branch_saliency[0]  = Float.NaN;

        if (verbose) appendFuzzification(ncc_2, likelihood_2, smoothness_2);

        branchStrength(ncc_2, likelihood_2, smoothness_2, t);
        float b2_is_off 	= t[0];
        float b2_is_none 	= t[1];
        float b2_is_on 		= t[2];
        branch_saliency[1]  = Float.NaN;

        if (verbose) appendFuzzification(ncc_3, likelihood_3, smoothness_3);

        branchStrength(ncc_3, likelihood_3, smoothness_3, t);
        float b3_is_off 	= t[0];
        float b3_is_none 	= t[1];
        float b3_is_on 		= t[2];
        branch_saliency[2]  = Float.NaN;

        if (verbose) appendFuzzification(ncc_4, likelihood_4, smoothness_4);

        branchStrength(ncc_4, likelihood_4, smoothness_4, t);
        float b4_is_off 	= t[0];
        float b4_is_none 	= t[1];
        float b4_is_on 		= t[2];
        branch_saliency[3]  = Float.NaN;

        Arrays.fill(agg2, 0);
        float mu;
        float[] cur;

        mu = min(b1_is_off, 	b2_is_off, 	b3_is_off, 	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0000
        mu = min(b1_is_off, 	b2_is_off, 	b3_is_off, 	b4_is_on);       cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0001 (end)
        mu = min(b1_is_off, 	b2_is_off, 	b3_is_on, 	b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0010 (end)
        mu = min(b1_is_off, 	b2_is_off, 	b3_is_on,  	b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0011

        mu = min(b1_is_off, 	b2_is_on,   b3_is_off,  b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 0100 (end)
        mu = min(b1_is_off, 	b2_is_on, 	b3_is_off,  b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0101
        mu = min(b1_is_off, 	b2_is_on, 	b3_is_on,  	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 0110
        mu = min(b1_is_off, 	b2_is_on, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 0111 (bif)

        mu = min(b1_is_on, 		b2_is_off, 	b3_is_off, 	b4_is_off);      cur = fi_endpoint(mu);  accumulate(cur, agg2); // 1000 (end)
        mu = min(b1_is_on,  	b2_is_off,  b3_is_off,  b4_is_on);       cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1001
        mu = min(b1_is_on,  	b2_is_off, 	b3_is_on,  	b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1010
        mu = min(b1_is_on,  	b2_is_off, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1011 (bif)

        mu = min(b1_is_on,  	b2_is_on, 	b3_is_off,  b4_is_off);      cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1100
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_off,  b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1101 (bif)
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_on,  	b4_is_off);      cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1110 (bif)
        mu = min(b1_is_on,  	b2_is_on, 	b3_is_on,  	b4_is_on);       cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1111 (bif)

//        mu = max(b1_is_none,    b2_is_none, b3_is_none, b4_is_none);     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        // don't allow 2 none, similarly as 1 none is not allowed when there are three
        mu = min(b1_is_none,    b2_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        mu = min(b1_is_none,    b3_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        mu = min(b1_is_none,    b4_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        mu = min(b2_is_none,    b3_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        mu = min(b2_is_none,    b4_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);
        mu = min(b3_is_none,    b4_is_none                        );     cur = fi_nonpoint(mu);  accumulate(cur, agg2);

//        mu = min(b1_is_none,	b2_is_on, 	b3_is_on, 	b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // x111 (bif)
//        mu = min(b1_is_on,		b2_is_none, b3_is_on, 	b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 1x11 (bif)
//        mu = min(b1_is_on,		b2_is_on, 	b3_is_none, b4_is_on);		 cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 11x1 (bif)
//        mu = min(b1_is_on,		b2_is_on, 	b3_is_on,   b4_is_none);     cur = fi_bifpoint(mu);  accumulate(cur, agg2); // 111x (bif)

//        mu = min(b1_is_none,	1-b2_is_on, 1-b3_is_on,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // x1_1_1_
//        mu = min(1-b1_is_on,	b2_is_none, 1-b3_is_on,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_x1_1_
//        mu = min(1-b1_is_on,	1-b2_is_on, b3_is_none,  1-b4_is_on);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_1_x1_
//        mu = min(1-b1_is_on,	1-b2_is_on, 1-b3_is_on,  b4_is_none);    cur = fi_nonpoint(mu);  accumulate(cur, agg2); // 1_1_1_x

        float defuzz = get_agg2_centroid();

        output_endpoint_nonpoint_bifpoint[0] = (float) Math.exp(-Math.pow(defuzz - output2_endpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[1] = (float) Math.exp(-Math.pow(defuzz - output2_nonpoint, 2) / (2 * Math.pow(output_sigma, 2)));
        output_endpoint_nonpoint_bifpoint[2] = (float) Math.exp(-Math.pow(defuzz - output2_bifpoint, 2) / (2 * Math.pow(output_sigma, 2)));

        // additional tweak - don't allow it to be higher than agg2 value at the defuzz point
        float limit = get_agg2(defuzz);

        output_endpoint_nonpoint_bifpoint[0] = Math.min(limit, output_endpoint_nonpoint_bifpoint[0]);
        output_endpoint_nonpoint_bifpoint[1] = Math.min(limit, output_endpoint_nonpoint_bifpoint[1]);
        output_endpoint_nonpoint_bifpoint[2] = Math.min(limit, output_endpoint_nonpoint_bifpoint[2]);


        String title =  "END=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[0],2) +
                "NON=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[1],2) +
                "JUN=" + IJ.d2s(output_endpoint_nonpoint_bifpoint[2],2);
        if (verbose) appendAgg2(defuzz, title);

    }

    private float get_agg1_centroid()
    {
        float centroid = 0;
        float norm = 0;
        for (int i=0; i<N1; i++) {
            centroid += x1[i]*agg1[i];
            norm += agg1[i];
        }

        if (norm>Float.MIN_VALUE) return centroid / norm;
        else return Stat.average(x1);

    }

    private float get_agg2_centroid()
    {
        float centroid = 0;
        float norm = 0;
        for (int i=0; i<N2; i++) {
            centroid += x2[i]*agg2[i];
            norm += agg2[i];
        }

        if (norm>Float.MIN_VALUE) return centroid / norm;
        else return Stat.average(x2);

    }

    private float get_agg2(float atX)
    {
        float step = (output2_range_ended-output2_range_start) / (N2-1);

        float out = Float.NaN;

        for (int i = 0; i < N2; i++) {
            float diff = Math.abs(atX-x2[i]);
            if (diff <= step/2f) {
                out = agg2[i];
                break;
            }
        }

        return out;

    }

    private static final float min(float in1, float in2)
    {
        return Math.min(in1, in2);
    }

    private static final float min(float in1, float in2, float in3)
    {
        return Math.min(in1, min(in2, in3));
    }

    private static final float min(float in1, float in2, float in3, float in4)
    {
        return Math.min(in1, min(in2, in3, in4));
    }

    private static final void accumulate(float[] values, float[] accumulator)
    {
        for (int i=0; i<values.length; i++) if (values[i]>accumulator[i]) accumulator[i] = values[i];
    }

    private static final float max(float in1, float in2)
    {
        return Math.max(in1, in2);
    }

    private static final float max(float in1, float in2, float in3)
    {
        return Math.max(in1, max(in2, in3));
    }

    private static final float max(float in1, float in2, float in3, float in4)
    {
        return Math.max(in1, max(in2, in3, in4));
    }

    private void appendAgg1(float defuzz_to_plot, String title)
    {
        Plot p = new Plot("", "", "");
        p.setLimits(x_min, x_max, 0 ,1);
        p.addPoints(x1, agg1, Plot.LINE);
        p.draw();
        p.setColor(Color.RED);
        p.setLineWidth(4);
        p.addPoints(new float[]{defuzz_to_plot, defuzz_to_plot}, new float[]{0, 1}, Plot.LINE);
        p.draw();
        fls_steps.addSlice(title, p.getProcessor());
    }

    private void appendAgg2(float defuzz_to_plot, String title)
    {
        Plot p = new Plot("", "", "");
        p.setLimits(x_min, x_max, 0 , 1);
        p.addPoints(x2, agg2, Plot.LINE);
        p.draw();
        p.setColor(Color.RED);
        p.setLineWidth(4);
        p.addPoints(new float[]{defuzz_to_plot, defuzz_to_plot}, new float[]{0, 1}, Plot.LINE);
        p.draw();
        fls_steps.addSlice(title, p.getProcessor());
    }

    private void appendFuzzification(float ncc_i, float likelihood_i, float smoothness_i)
    {
        int NN = 51;
        float[] xx = new float[NN];
        float[] yy_high = new float[NN]; // all the inputs will be fuzzified as high/low membership degree
        float[] yy_low = new float[NN];

        // input 1: ncc
        for (int ii=0; ii<xx.length; ii++) xx[ii] = ncc_start + ii * ( (ncc_end-ncc_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = (h_ncc_high(xx[ii])<=h_ncc_high(ncc_i))? h_ncc_high(xx[ii]) : h_ncc_high(ncc_i) ;
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = (h_ncc_low(xx[ii]) <=h_ncc_low(ncc_i))? h_ncc_low(xx[ii])  : h_ncc_low(ncc_i);
        Plot p = new Plot("", "", "");
        p.setLimits(ncc_start, ncc_end, 0, 1);
        p.setColor(Color.RED);
        p.addPoints(xx, yy_high, Plot.LINE);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        p.draw();
        p.setColor(Color.DARK_GRAY);
        p.setLineWidth(4);
        p.addPoints(new float[]{ncc_i, ncc_i}, new float[]{0, 1}, Plot.LINE);
        p.draw();
        String tit_ncc = "NCC="+IJ.d2s(ncc_i,2)+" (LOW=" + IJ.d2s(h_ncc_low(ncc_i),1) + ",HIGH=" + IJ.d2s(h_ncc_high(ncc_i),1) + ")";
        fls_steps.addSlice(tit_ncc, p.getProcessor());

        // input 2: likelihood
        for (int ii=0; ii<xx.length; ii++) xx[ii] = likelihood_start + ii * ( (likelihood_end-likelihood_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = (h_likelihood_high(xx[ii])<=h_likelihood_high(likelihood_i))? h_likelihood_high(xx[ii]) : h_likelihood_high(likelihood_i);
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = (h_likelihood_low(xx[ii])<=h_likelihood_low(likelihood_i))? h_likelihood_low(xx[ii]) : h_likelihood_low(likelihood_i);
        p = new Plot("", "", "");
        p.setLimits(likelihood_start, likelihood_end, 0, 1);
        p.setColor(Color.RED);
        p.addPoints(xx, yy_high, Plot.LINE);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        p.draw();
        p.setColor(Color.DARK_GRAY);
        p.setLineWidth(4);
        p.addPoints(new float[]{likelihood_i, likelihood_i}, new float[]{0, 1}, Plot.LINE);
        p.draw();
        String tit_lhood = "LHOOD="+IJ.d2s(likelihood_i,2)+" (LOW=" + IJ.d2s(h_likelihood_low(likelihood_i),1) + ",HIGH=" + IJ.d2s(h_likelihood_high(likelihood_i),1) + ")";
        fls_steps.addSlice(tit_lhood, p.getProcessor());

        // input 3: smoothness
        for (int ii=0; ii<xx.length; ii++) xx[ii] = smoothness_start + ii * ( (smoothness_end-smoothness_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = (h_smoothness_high(xx[ii])<=h_smoothness_high(smoothness_i))? h_smoothness_high(xx[ii]) : h_smoothness_high(smoothness_i);
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = (h_smoothness_low(xx[ii])<=h_smoothness_low(smoothness_i))? h_smoothness_low(xx[ii]) : h_smoothness_low(smoothness_i);
        p = new Plot("", "", "");
        p.setLimits(smoothness_start, smoothness_end, 0, 1);
        p.setColor(Color.RED);
        p.addPoints(xx, yy_high, Plot.LINE);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        p.draw();
        p.setColor(Color.DARK_GRAY);
        p.setLineWidth(4);
        p.addPoints(new float[]{smoothness_i, smoothness_i}, new float[]{0, 1}, Plot.LINE);
        p.draw();
        String tit_sthness = "SMTH="+IJ.d2s(smoothness_i,2)+",SMTH_LOW/HIGH"+smoothness_low+""+smoothness_high+", (LOW=" + IJ.d2s(h_smoothness_low(smoothness_i),1) + ",HIGH=" + IJ.d2s(h_smoothness_high(smoothness_i),1) + ")";
        fls_steps.addSlice(tit_sthness, p.getProcessor());

    }

    public void showFuzzification()
    {

        ImageStack is_out = new ImageStack(528, 255);
        int NN = 51; // number of points for plotting

        // ncc input fuzzification
        float[] xx = new float[NN];
        float[] yy_high = new float[NN]; // all the inputs will be fuzzified as high/low membership degree
        float[] yy_low = new float[NN];

        // input 1: ncc
        for (int ii=0; ii<xx.length; ii++) xx[ii] = ncc_start + ii * ( (ncc_end-ncc_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = h_ncc_high(xx[ii]);
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = h_ncc_low(xx[ii]);
        Plot p = new Plot("ncc", "", "", xx, yy_high, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        is_out.addSlice("NCC", p.getProcessor());

        // input 2: likelihood
        for (int ii=0; ii<xx.length; ii++) xx[ii] = likelihood_start + ii * ( (likelihood_end-likelihood_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = h_likelihood_high(xx[ii]);
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = h_likelihood_low(xx[ii]);

        p = new Plot("likelihood", "", "", xx, yy_high, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        is_out.addSlice("LIKELIHOOD", p.getProcessor());

        // input 3: smoothness
        for (int ii=0; ii<xx.length; ii++) xx[ii] = smoothness_start + ii * ( (smoothness_end-smoothness_start) / (NN-1));
        for (int ii=0; ii<xx.length; ii++) yy_high[ii] = h_smoothness_high(xx[ii]);
        for (int ii=0; ii<xx.length; ii++) yy_low[ii] = h_smoothness_low(xx[ii]);

        p = new Plot("smoothness", "", "", xx, yy_high, Plot.LINE);
        p.setColor(Color.RED);
        p.draw();
        p.setColor(Color.BLUE);
        p.addPoints(xx, yy_low, Plot.LINE);
        is_out.addSlice("SMOOTHNESS", p.getProcessor());

        new ImagePlus("FUZZIFICATION", is_out).show();

    }

    public void showDefuzzification()
    {
        ImageStack is_out = new ImageStack(528,255);

        Plot p = new Plot("", "", ""); p.setLimits(x_min, x_max, 0, 1);
        p.setColor(Color.RED);
        p.addPoints(x1, q_on, Plot.LINE);
        p.setColor(Color.GREEN);
        p.addPoints(x1, q_none, Plot.LINE);
        p.setColor(Color.BLUE);
        p.addPoints(x1, q_off, Plot.LINE);
        is_out.addSlice("OUT1", p.getProcessor());

        p = new Plot("", "", ""); p.setLimits(x_min, x_max, 0, 1);
        p.setColor(Color.RED);
        p.addPoints(x2, q_bifpoint, Plot.LINE);
        p.setColor(Color.GREEN);
        p.addPoints(x2, q_nonpoint, Plot.LINE);
        p.setColor(Color.BLUE);
        p.addPoints(x2, q_endpoint, Plot.LINE);
        is_out.addSlice("OUT2", p.getProcessor());

        new ImagePlus("DEFUZZIFICATION", is_out).show();

    }

}
