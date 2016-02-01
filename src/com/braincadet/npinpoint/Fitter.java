package com.braincadet.npinpoint;

import ij.IJ;
import ij.ImageStack;

import java.util.ArrayList;

/** fitting set of profile data to the template array */
public class Fitter {

    int patch_radial_len;
    int patch_tangen_len;

    // template profiles
    public ArrayList<float[]> 	templates; 						// this is 2d template actually
    public ArrayList<Float> 	templates_mean;
    public ArrayList<float[]>	templates_mean_subtr;
    public ArrayList<Float> 	templates_mean_subtr_sumsqr;
    public ArrayList<String>	template_legend;

    // middle
    float middle_idx,
            start_sigma, d_sigma, end_sigma, // slope
            start_width, d_width, end_width; // basic width

    public Fitter(int _radial_len, int _tangen_len, float _cross_sigma_ratio, boolean verbose) // _tangen_len has to be 2N+1 length
    {
        patch_radial_len = _radial_len;
        patch_tangen_len = _tangen_len;

        middle_idx = (patch_tangen_len-1) / 2f;
        // sigma
        start_sigma = _cross_sigma_ratio * patch_tangen_len * 0.7f;//*0.02f;
        d_sigma     = _cross_sigma_ratio * patch_tangen_len * 0.3f; // dummy
        end_sigma   = _cross_sigma_ratio * patch_tangen_len * 1.31f; // only one value of sigma
        // width
        start_width = 0;
        d_width     = patch_tangen_len;
        end_width   = 0;

        // variables will be stored in:
        // templates                                  	(vector)
        // templates_mean (used for NCC calculation)  	(scalar)
        // templates - templates_mean                 	(vector)
        // sum(templates - templates_mean)				(scalar)
        // templates_legend 							(string)
        templates 			= new ArrayList<float[]>();
        templates_mean 		= new ArrayList<Float>();
        templates_mean_subtr= new ArrayList<float[]>();
        templates_mean_subtr_sumsqr = new ArrayList<Float>();
        template_legend = new ArrayList<String>();

        // create templates
        for(float width = start_width; width<=end_width; width+=d_width) {

            for (float sigma = start_sigma; sigma <= end_sigma; sigma+=d_sigma) {

                float[] templates_element = new float[patch_radial_len*patch_tangen_len];

                float boundary_1 = Math.round(middle_idx - width/2);
                float boundary_2 = Math.round(middle_idx + width/2);

                for (int j = 0; j < patch_radial_len; j++) {
                    for (int i = 0; i < patch_tangen_len; i++) {

                        float val, d;
                        if (i < boundary_1) {
                            d = boundary_1 - i;
                            val = (float) Math.exp(-Math.pow(d, 2)/(2*sigma*sigma));
                        }
                        else if (i >= boundary_1 && i < boundary_2) {
                            val = 1;
                        }
                        else {
                            d = i - boundary_2;
                            val = (float) Math.exp(-Math.pow(d, 2)/(2*sigma*sigma));
                        }

                        templates_element[j*patch_tangen_len+i] = val;

                    }
                }

                templates.add(templates_element.clone());

                // check boundary elements are low enough
                //if (templates_element[0] < low_boundary && templates_element[vector_len-1] < low_boundary) {

                // calculate mean
                float mn = Stat.average(templates_element);
                templates_mean.add(mn);

                // subtract and store in the same array, add the sum as well
                float sum = 0;
                for (int aa = 0; aa < templates_element.length; aa++) {
                    templates_element[aa] = templates_element[aa] - mn;
                    sum += templates_element[aa] * templates_element[aa];
                }

                templates_mean_subtr.add(templates_element.clone());
                templates_mean_subtr_sumsqr.add(sum);
                template_legend.add("wdt="+ IJ.d2s(width, 2)+","+"sig="+IJ.d2s(sigma, 2));

            }

        }

        if (verbose) {
            System.out.println("created " + templates.size() + " template profiles.");
        }

    }

    public ImageStack getTemplates() {

        ImageStack is_out = new ImageStack(patch_tangen_len, patch_radial_len);

        for (int aa= 0; aa<templates.size(); aa++) { // plot each template
            is_out.addSlice(template_legend.get(aa), templates.get(aa));
        }

        return is_out;

    }

    public float[] getTemplate(int template_index) {

        return templates.get(template_index);

    }

    public float[] fit(float[] profile) {

        // returns the fitting result: (index of the profile, fitting score)
        float[] out = new float[2];
        out[0] = Float.NaN;
        out[1] = Float.NEGATIVE_INFINITY; // looking for  highest ncc

        // loop the templates
        for (int i=0; i<templates.size(); i++) {

            // calculate score
            float curr_score = ncc(profile, templates_mean_subtr.get(i), templates_mean_subtr_sumsqr.get(i));
            if (curr_score > out[1]) {out[0]=i; out[1]=curr_score;} // first is index, second is ncc score

        }

        return out;

    }

    private float ncc(float[] f, float[] t_tM, float sumsqr_t_tM)  // f needs to have same length as template
    {

        // mean(patch)
        float f_mean = f[0];
        for (int i=1; i < f.length; i++) f_mean += f[i];
        f_mean = f_mean / (float)f.length;

        // sumsqr(patch-mean(patch))
        float sc = 0;
        float f_sub_f_mean_sumsqr = 0;

        for (int aa=0; aa<(patch_radial_len*patch_tangen_len); aa++) {
            sc += (f[aa]-f_mean) * t_tM[aa]; // important that input f is vector_len length
            f_sub_f_mean_sumsqr += (f[aa]-f_mean) * (f[aa]-f_mean);
        }

        float ncc_val;

        if (f_sub_f_mean_sumsqr>Float.MIN_VALUE) {
            ncc_val = sc / (float) Math.sqrt(f_sub_f_mean_sumsqr * sumsqr_t_tM);
        }
        else {
            ncc_val = 0;
        }
        return ncc_val;

    }

}
