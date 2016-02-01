package com.braincadet.npinpoint;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Line;
import ij.gui.OvalRoi;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;

public class Detector {

    String 		image_dir;
    String		image_name;
    float[][] 	inimg_xy;               			// store input image as an array

    // parameters
    float       	s;
    float 			sigma_ratio; 					// sigma = sigma_ratio * D
    float[]         D;
    float       	ncc_high;
    float       	ncc_low;

    float 			likelihood_high;
    float 			likelihood_low;

    float           smoothness_high;            	// smoothness_high is actually lower than smoothness_low
    float           smoothness_low;             	// they are automatically calculated fro mthe distribution of smoothness values

    float[]			dsens;                          // detection sensitivities, score thresholds read from csv argument string

    float			output_sigma = 0.45f;

    float 			output_membership_th;       	// based on k and output_sigma

    int         	M = 1;
    float       	minCos = -.5f;
    float			k = 0.1f;                   	// hardcoded, defines the strictness of the threshold for out membership function, k higher => means more strictness

    // hardcoded (for grid parameter testing - does not make sense to examine more)
    int 			maxNrRegions = 500;            // max nr regions for CP type (more than that is discarded)
    float 			ratioOnPix	= 0.05f;		   // pixels that were belonging to CP regions cannot be more than ratio of the total area

    private boolean save_midresults = false;
    String          midresults_dir = "";
    public boolean  auto_smoothness = false;

    public boolean  do_junctions  = true;           // control flags
    public boolean  do_endpoints  = true;
    public boolean  do_directions = false;

    Fuzzy sample_fls = null; 				        // will be used for user interaction (to simulate detection for features extracted at single locations)

    EntrTh et = new EntrTh();                       // will be used to threshold the scores after regularisation

    float[] kernel = new float[9];     				// for score regularisation
    public int Nreg = 2; 							// how many times it is done
    public float ang_deg = 30f;  					// ms parameter

    int         CPU_NR;

    // OUTPUT
    ImagePlus ip_exporter = new ImagePlus(); 		// to wrap up the processors for saving & exporting

    FloatProcessor map_scores_end = null;			// scores map used to delineate the regions
    FloatProcessor map_scores_jun = null;           // store 2d grid with detection scores (size of the image itself)

    // Note: thresholding the scores with set of threshold levels (csv argument) results in list of following items (one per threshold):
    // map_region_end/jun, cumm_regions_end/jun, cumm_directions_end/jun, detected_regions
    // base list will be per detection sensitivity
    ArrayList<String>        output_dir_name; 		    // parameter coded output folder name

    ArrayList<ByteProcessor> map_region_end = null;     // region binary maps that update throughout the scales
    ArrayList<ByteProcessor> map_region_jun = null;     // pixel at (x,y) assigned to one region

    ArrayList<ByteProcessor> cumm_regions_end = null;   // after accumulating booleans throughout several scales
    ArrayList<ByteProcessor> cumm_regions_jun = null;   // still binary maps

    ArrayList<ArrayList[][]> cumm_directions_end = null;       // after accumulating the directions throughout several scales
    ArrayList<ArrayList[][]> cumm_directions_jun = null;       // accumulates 2d directions for each location

    ArrayList<ArrayList<CritpointRegion>> detected_regions; 	// list with the detections

    private static int plotw = new Plot("","","", new float[1], new float[1]).getSize().width;
    private static int ploth = new Plot("","","", new float[1], new float[1]).getSize().height;

    public Detector(
            ImagePlus 	ip_load,
            float 		_s,
            float[] 	_D,
            float 		_sigma_ratio,
            float 		_ncc_high,
            float 		_ncc_low,
            float 		_likelihood_high,
            float 		_likelihood_low,
            float      _smoothness_high,
            float      _smoothness_low,
            float[] 	_dsens
    )
    {

        image_dir = ip_load.getOriginalFileInfo().directory; 			//  + File.separator  + image_name
        image_name = ip_load.getShortTitle();

        inimg_xy = new float[ip_load.getWidth()][ip_load.getHeight()]; 	// x~column, y~row
        if (ip_load.getType()== ImagePlus.GRAY8) {
            byte[] read = (byte[]) ip_load.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = (float) (read[idx] & 0xff);
            }
        }
        else if (ip_load.getType()==ImagePlus.GRAY32) {
            float[] read = (float[]) ip_load.getProcessor().getPixels();
            for (int idx=0; idx<read.length; idx++) {
                inimg_xy[idx%ip_load.getWidth()][idx/ip_load.getWidth()] = read[idx];
            }
        }
        else {
            IJ.log("image type not recognized");
            return;
        }

        // detection scale define
        s = _s;
        D = _D;
        sigma_ratio = _sigma_ratio;

        // detection parameters
        ncc_high = _ncc_high;
        ncc_low = _ncc_low;
        likelihood_high = _likelihood_high;
        likelihood_low = _likelihood_low;
        smoothness_high = _smoothness_high;
        smoothness_low = _smoothness_low;

        dsens = _dsens;
        if (dsens.length==0) {
            IJ.log("no detection sensitivity thresholds provided");
            return;
        }

        // sort ascending
        Arrays.sort(dsens);

        String		Dlist = ""; // for the out directory
        for (int i=0; i<_D.length; i++) Dlist += IJ.d2s(_D[i], 1) + ((i == _D.length - 1) ? "" : ","); // take care that it is a valid folder name signal ',._' are valid

        output_membership_th = (float) Math.exp(-(0.5f*0.5f)/(2*output_sigma*output_sigma)); // 0.5 is middle of the output range
        output_membership_th = (float) (1 - Math.pow(output_membership_th,k) * (1-output_membership_th)); // h1 = 1 - h^k * (1-h)

        output_dir_name = new ArrayList<String>(dsens.length);
        for (int i = 0; i < dsens.length; i++) {
            String pattern = image_dir+String.format(
                    "DET.D.nncH.nccL.lhoodH.lhoodL.sthH.sthL.dsens_"+
                            "%s_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f_%.2f",
                    Dlist,
                    ncc_high,
                    ncc_low,
                    likelihood_high,
                    likelihood_low,
                    smoothness_high,
                    smoothness_low,
                    dsens[i]
            );
            output_dir_name.add(i, pattern);

        }

        midresults_dir = image_dir+image_name + "_midresults" + File.separator;

        CPU_NR = Runtime.getRuntime().availableProcessors() + 1; //5;

        map_scores_end = new FloatProcessor(ip_load.getWidth(), ip_load.getHeight());
        map_scores_jun = new FloatProcessor(ip_load.getWidth(), ip_load.getHeight());

        // allocate outputs
        map_region_end = new ArrayList<ByteProcessor>(dsens.length);
        map_region_jun = new ArrayList<ByteProcessor>(dsens.length);
        // no need to allocate these further, just add() once they are rteturned

        // initialize, these are still boolean matrices, but need to be allocated
        cumm_regions_end = new ArrayList<ByteProcessor>(dsens.length);
        cumm_regions_jun = new ArrayList<ByteProcessor>(dsens.length);
        for (int i = 0; i < dsens.length; i++) {
            cumm_regions_end.add(i, new ByteProcessor(ip_load.getWidth(), ip_load.getHeight()));
            cumm_regions_jun.add(i, new ByteProcessor(ip_load.getWidth(), ip_load.getHeight()));
        }

        cumm_directions_end = new ArrayList<ArrayList[][]>(dsens.length);
        cumm_directions_jun = new ArrayList<ArrayList[][]>(dsens.length);
        for (int i = 0; i < dsens.length; i++) {

            ArrayList[][] tmp_end = new ArrayList[ip_load.getWidth()][ip_load.getHeight()];
            ArrayList[][] tmp_jun = new ArrayList[ip_load.getWidth()][ip_load.getHeight()];

            for (int ww = 0; ww < ip_load.getWidth(); ww++) {
                for (int hh = 0; hh < ip_load.getHeight(); hh++) {
                    tmp_end[ww][hh] = null;
                    tmp_jun[ww][hh] = null;
                }
            }

            cumm_directions_end.add(i, tmp_end);
            cumm_directions_jun.add(i, tmp_jun);

        }

        Arrays.fill(kernel, 1/9f);

        detected_regions = new ArrayList<ArrayList<CritpointRegion>>(dsens.length);
        for (int i = 0; i < dsens.length; i++) detected_regions.add(i, new ArrayList<CritpointRegion>());

        // used for visualization, normally features are extracted using Fuzzy instance in FuzzyDetector
        sample_fls = new Fuzzy( // default initialization with given smoothness (yet to calculate auto)
                ncc_high, // ncc
                ncc_low,
                likelihood_high, // lhood
                likelihood_low,
                smoothness_high, // smoothness default
                smoothness_low,
                output_sigma  // std output membership functions - defines separation
        );

    }

    public void saveMidresults(boolean flag) {
        // will put the flag on and create the output folder for midresults

        save_midresults = flag;

        if (save_midresults) {  // if set to true create the folder as well
            File f1 = new File(midresults_dir);
            if (!f1.exists()) {
                f1.mkdirs();
            }
            else { // clean it, reset with each exec
                File[] files = f1.listFiles();
                for (File file : files)
                {
                    if (!file.delete())
                        IJ.log("Failed to delete " + file);
                }
            }
        }


    }

    public void run()
    {

        // create output directories
        for (int i = 0; i < output_dir_name.size(); i++) {
            File f = new File(output_dir_name.get(i));
            if (!f.exists()) f.mkdirs();
        }

        long t1, t2;
        t1 = System.currentTimeMillis();

        for (int didx = 0; didx<D.length; didx++) { // loop scales D[]

            Sphere2D sph2d = new Sphere2D(D[didx], s, sigma_ratio);
            if (save_midresults) {
                IJ.saveAs(sph2d.showSampling(), "Tiff", midresults_dir+"sampling_"+D[didx]+".tif");
                IJ.saveAs(sph2d.showWeights(),  "Tiff", midresults_dir+"weights_"+D[didx]+".tif");
            }
            /********************************************************************/
            IJ.log("Masker...");
            float new_masker_radius = 1.5f*sph2d.getOuterRadius();   	// important that it is outer radius of the sphere
            float new_masker_percentile = 50;                   		// used to have these two as argument but not necessary
            Masker.loadTemplate(
                    inimg_xy,
                    (int)Math.ceil(new_masker_radius),
                    new_masker_radius,
                    new_masker_percentile); //image, margin, check, percentile
            int totalLocs = inimg_xy.length * inimg_xy[0].length;
            Masker ms_jobs[] = new Masker[CPU_NR];
            for (int i = 0; i < ms_jobs.length; i++) {
                ms_jobs[i] = new Masker(i*totalLocs/CPU_NR,  (i+1)*totalLocs/CPU_NR);
                ms_jobs[i].start();
            }
            for (int i = 0; i < ms_jobs.length; i++) {
                try {
                    ms_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            Masker.defineThreshold();
            Masker.formRemainingOutputs();
            if (save_midresults) {
                ImagePlus mask = new ImagePlus("mask", Masker.getMask());
                IJ.saveAs(mask, "Tiff", midresults_dir+"mask_"+D[didx]+".tif");
            }
            /********************************************************************/
            IJ.log("Profiler...");
            Profiler.loadTemplate(
                    sph2d,
                    Masker.i2xy,
                    Masker.xy2i,
                    inimg_xy);
            int totalProfileComponents = sph2d.getProfileLength();
            Profiler pf_jobs[] = new Profiler[CPU_NR];
            for (int i = 0; i < pf_jobs.length; i++) {
                pf_jobs[i] = new Profiler(i*totalProfileComponents/CPU_NR,  (i+1)*totalProfileComponents/CPU_NR);
                pf_jobs[i].start();
            }
            for (int i = 0; i < pf_jobs.length; i++) {
                try {
                    pf_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            /********************************************************************/
            IJ.log("ProfileSpread...");
            ProfileSpread.loadTemplate(
                    Masker.i2xy,
                    Masker.xy2i,
                    Profiler.prof2,
                    inimg_xy.length,
                    inimg_xy[0].length);
            int totalProfileLocations = Profiler.prof2.length;
            ProfileSpread pv_jobs[] = new ProfileSpread[CPU_NR];
            for (int i = 0; i < pv_jobs.length; i++) {
                pv_jobs[i] = new ProfileSpread(i*totalProfileLocations/CPU_NR,  (i+1)*totalProfileLocations/CPU_NR);
                pv_jobs[i].start();
            }
            for (int i = 0; i < pv_jobs.length; i++) {
                try {
                    pv_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            ProfileSpread.threshold();
            if (save_midresults) {
                ImagePlus testip = new ImagePlus("", ProfileSpread.getMask());
                IJ.saveAs(testip, "Tiff", midresults_dir+"mask_profile_"+D[didx]+".tif");
            }
            IJ.log(" " + IJ.d2s((ProfileSpread.getNrCritpointCandidates() * 100f) / (inimg_xy.length * inimg_xy[0].length), 2) + " % candidates...");
            /********************************************************************/
            IJ.log("PeakExtractor...");
            PeakExtractor.loadTemplate(sph2d, Masker.i2xy, Profiler.prof2, inimg_xy, Masker.xy2i);
            int totalPeakExtrComponents = Profiler.prof2.length; // number of profiles == number of locations i2xy.length
            PeakExtractor pe_jobs[] = new PeakExtractor[CPU_NR];
            for (int i = 0; i < pe_jobs.length; i++) {
                pe_jobs[i] = new PeakExtractor(i*totalPeakExtrComponents/CPU_NR, (i+1)*totalPeakExtrComponents/CPU_NR);
                pe_jobs[i].start();
            }
            for (int i = 0; i < pe_jobs.length; i++) {
                try {
                    pe_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            //PeakExtractor.getCircStat().show();
            /********************************************************************/
            IJ.log("Delineator...");
            Delineator.loadTemplate(
                    Masker.i2xy,
                    Masker.xy2i,
                    ProfileSpread.profile_diverse,
                    PeakExtractor.peaks_i,
                    PeakExtractor.peaks_w,
                    inimg_xy,
                    sph2d,	//D[didx],
                    M,
                    minCos
            );
            int totalDelineationComponents = PeakExtractor.peaks_i.length;
            Delineator dl_jobs[] = new Delineator[CPU_NR];
            for (int i = 0; i < dl_jobs.length; i++) {
                dl_jobs[i] = new Delineator(i*totalDelineationComponents/CPU_NR, (i+1)*totalDelineationComponents/CPU_NR);
                dl_jobs[i].start();
            }
            for (int i = 0; i < dl_jobs.length; i++) {
                try {
                    dl_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            if (save_midresults) {  // TODO
                //new ImagePlus("", Delineator.getSmoothnessDistribution(64)).show();
            }

            if (auto_smoothness) {
                int percentile = 90;
                float sensitivity = 0.5f;
                smoothness_high = Delineator.getSmoothnessPercentile(percentile);
                smoothness_low = sensitivity * smoothness_high;
                sample_fls = new Fuzzy(// redefine flas used later for simulating, for visualizations (should be the same as the one in FuzzyDetector.run())
                        ncc_high, // ncc
                        ncc_low,
                        likelihood_high, // lhood
                        likelihood_low,
                        smoothness_high, // smoothness redefined
                        smoothness_low,
                        output_sigma  // std output membership functions - defines separation
                );
            }
            /********************************************************************/
            IJ.log("Ncc...");
            Ncc.loadTemplate(
                    Delineator.xy2,
                    Delineator.vxy2,
                    Delineator.L,
                    Delineator.dim,
                    Delineator.samplingStep,
                    inimg_xy,
                    sigma_ratio
            );
            int totalNccExtractionComponents = Delineator.xy2.length;
            Ncc ncc_jobs[] = new Ncc[CPU_NR];
            for (int i = 0; i < ncc_jobs.length; i++) {
                ncc_jobs[i] = new Ncc(i*totalNccExtractionComponents/CPU_NR, (i+1)*totalNccExtractionComponents/CPU_NR);
                ncc_jobs[i].start();
            }
            for (int i = 0; i < ncc_jobs.length; i++) {
                try {
                    ncc_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
            if (save_midresults) {
                ImagePlus testip = new ImagePlus("", Ncc.getTemplates());
                IJ.saveAs(testip, "Tiff", midresults_dir+"fitter_templates_"+D[didx]+".tif");
            }

            /********************************************************************/
            IJ.log("FuzzyDetector...");
            FuzzyDetector.loadTemplate(
                    Ncc.scores,
                    PeakExtractor.peaks_lhood,
                    Delineator.smoothness,

                    ncc_high,        // user
                    ncc_low,         // user

                    likelihood_high,  // user
                    likelihood_low,   // user

                    smoothness_high,  // automatic or manual
                    smoothness_low,   // automatic or manual

                    output_sigma
            );
            int totalFuzzyDetectorComponents = PeakExtractor.peaks_lhood.length;
            FuzzyDetector fd_jobs[] = new FuzzyDetector[CPU_NR];
            for (int i = 0; i < fd_jobs.length; i++) {
                fd_jobs[i] = new FuzzyDetector(i*totalFuzzyDetectorComponents/CPU_NR, (i+1)*totalFuzzyDetectorComponents/CPU_NR);
                fd_jobs[i].start();
            }
            for (int i = 0; i < fd_jobs.length; i++) {
                try {
                    fd_jobs[i].join();
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }

            /********************************************************************/
            if (do_endpoints) {

                IJ.log("AppendEndpoints...");
                fillUp(map_scores_end, 0); // reset before each scale
                fillUp(map_scores_end, FuzzyDetector.endpoint_score, Masker.i2xy); // 2d

                if (save_midresults) {
                    ip_exporter.setProcessor(map_scores_end);
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_scores_end_"+D[didx]+".tif");
                }

                // regularization
                for (int i = 0; i < Nreg; i++) map_scores_end.convolve(kernel, 3, 3); // regularization

                if (save_midresults) {
                    ip_exporter.setProcessor(map_scores_end);
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_scores_end_reg_"+D[didx]+".tif");
                }

                map_region_end.clear();

                for (int i = 0; i < dsens.length; i++) {  // for all dsens thresholds

                    map_region_end.add(i, et.runMaxEntropyThreshold(map_scores_end, dsens[i])); // map_scores_end -> map_region_end, use maximum entropy to threshold ~1

                    if (save_midresults) {
                        ip_exporter.setProcessor(map_region_end.get(i));
                        IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_region_end_D="+D[didx] + ",dsens=" + IJ.d2s(dsens[i],2) + ".tif");
                    }

                    // append  to the all-scale region map (using "OR")
                    appendLogicalOr(cumm_regions_end.get(i), map_region_end.get(i));

                    if (do_directions) { // append all-scale direction map
                        appendDirections(cumm_directions_end.get(i), map_region_end.get(i),
                                Profiler.xy2i, PeakExtractor.i2xy, PeakExtractor.peaks_i, FuzzyDetector.branch_score, output_membership_th);
                    }// otherwise it stays with all the components initialized to null <[][]=null>, use cumm_directions_end<> only with do_directions label

                }

            }
            /********************************************************************/
            if (do_junctions) {

                IJ.log("AppendJunctions...");
                fillUp(map_scores_jun, 0); // reset before each scale
                fillUp(map_scores_jun, FuzzyDetector.junction_score, Masker.i2xy); // 2d

                if (save_midresults) {
                    ip_exporter.setProcessor(map_scores_jun);
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_scores_jun_"+D[didx]+".tif");
                }

                // regularization
                for (int i = 0; i < Nreg; i++) map_scores_jun.convolve(kernel, 3, 3); // regularization

                if (save_midresults) {
                    ip_exporter.setProcessor(map_scores_jun);
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_scores_jun_reg_"+D[didx]+".tif");
                }

                map_region_jun.clear();

                for (int i = 0; i < dsens.length; i++) {  // for all thresholds

                    map_region_jun.add(i, et.runMaxEntropyThreshold(map_scores_jun, dsens[i])); // map_scores_end -> map_region_end, use maximum entropy to threshold ~1

                    if (save_midresults) {
                        ip_exporter.setProcessor(map_region_jun.get(i));
                        IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"map_region_jun_D="+D[didx] + ",dsens=" + IJ.d2s(dsens[i],2) + ".tif");
                    }

                    // append  to the all-scale region map (using "OR")
                    appendLogicalOr(cumm_regions_jun.get(i), map_region_jun.get(i));


                    if (do_directions) { // append all-scale direction map
                        appendDirections(cumm_directions_jun.get(i), map_region_jun.get(i),
                                Profiler.xy2i, PeakExtractor.i2xy, PeakExtractor.peaks_i, FuzzyDetector.branch_score, output_membership_th);
                    }// otherwise it stays with all the components initialized to null <[][]=null>, use cumm_directions_jun<> only with do_directions label

                }

            }

            //IJ.log(" " + didx + "/" + (D.length-1) );

        } // loop D[]

        t2 = System.currentTimeMillis();
        IJ.log("\nfinished. " + ((t2 - t1) / 1000f) + "sec.");
        t1 = System.currentTimeMillis();

        /** addition */
//        IJ.log("testing... " + cumm_regions_end.size() + " || " + cumm_regions_jun.size());
//        for (int i = 0; i < dsens.length; i++) IJ.log(i + " : " + dsens[i] + " , ");

        // dsens is ascending!!!! redefine the cumm_regions_*
//        for (int i = dsens.length-2; i >= 0; i--) {
        for (int i = 1; i < dsens.length; i++) {
            appendLogicalOr(cumm_regions_end.get(i), cumm_regions_end.get(i-1));
            appendLogicalOr(cumm_regions_jun.get(i), cumm_regions_jun.get(i-1));
        }

        if (do_endpoints) {

            IJ.log("EndpointRegions... ");

            for (int i = 0; i <cumm_regions_end.size(); i++) {

//                IJ.log("detected_regions.get("+i+")");
//                IJ.log("before - " + detected_regions.get(i).size());

                if (save_midresults) {
                    ip_exporter.setProcessor(cumm_regions_end.get(i));
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir + "cumm_regions_end,dsens=" + IJ.d2s(dsens[i],2) + ".tif");
                }

                if (do_directions) {
                    IJ.log("(+MS to extract directions)... ");
                    appendDetectedRegions(
                            cumm_regions_end.get(i),
                            cumm_directions_end.get(i),
                            map_scores_end, CritpointRegion.RegionType.END, ang_deg,
                            detected_regions.get(i)); // append in detected_regions
                }
                else {
                    appendDetectedRegions( // no direction analysis
                            cumm_regions_end.get(i),
                            map_scores_end, CritpointRegion.RegionType.END,
                            detected_regions.get(i)
                    );
                }

//                IJ.log("after - " + detected_regions.get(i).size());

            }

        }

        if (do_junctions) {

            IJ.log("JunctionRegions... ");

            for (int i = 0; i <cumm_regions_jun.size(); i++) {

                if (save_midresults) {
                    ip_exporter.setProcessor(cumm_regions_jun.get(i));
                    IJ.saveAs(ip_exporter, "Tiff", midresults_dir+"cumm_regions_jun,dsens="+IJ.d2s(dsens[i],2)+".tif");
                }

                if (do_directions) {
                    IJ.log("(+MS to extract directions)... ");
                    appendDetectedRegions(
                            cumm_regions_jun.get(i),
                            cumm_directions_jun.get(i),
                            map_scores_jun, CritpointRegion.RegionType.BIF_CROSS, ang_deg,
                            detected_regions.get(i)); // append in detected_regions
                }
                else {
                    appendDetectedRegions( // no direction analysis
                            cumm_regions_jun.get(i),
                            map_scores_jun, CritpointRegion.RegionType.BIF_CROSS,
                            detected_regions.get(i)
                    );
                }


            }

        }

        t2 = System.currentTimeMillis();
        IJ.log("\ndone. " + ((t2 - t1) / 1000f) + "sec.");

    } // run()

    private static void fillUp(FloatProcessor critpoint2d, float[] critpoint1d, int[][] _i2xy) {
        for (int ll=0; ll<critpoint1d.length; ll++) {

            int x = _i2xy[ll][0];
            int y = _i2xy[ll][1];

            critpoint2d.setf(x, y, critpoint1d[ll]); //critpoint2d[x][y] = critpoint1d[ll];

        }
    }

    private static void fillUp(FloatProcessor template, float val)
    {
        float[] vals = (float[]) template.getPixels();
        for (int i = 0; i < vals.length; i++) {
            vals[i] = val;
        }
    }

    private static float getMaximum(FloatProcessor template)
    {

        float[] vals = (float[]) template.getPixels();
        float max = Float.NEGATIVE_INFINITY;
        for (int i = 0; i < vals.length; i++) {
            if (vals[i] > max) {
                max = vals[i];
            }
        }
        return max;

    }

    private static void threshold(FloatProcessor in, float th, ByteProcessor out)
    {
        float[] vals = (float[]) in.getPixels();
        for (int i = 0; i < vals.length; i++) {
            if (vals[i] >=th) {
                out.set(i, 255);
            }
            else {
                out.set(i, 0);
            }
        }

    }

    private static void appendLogicalOr(ByteProcessor base, ByteProcessor _to_append)
    {
        for (int i = 0; i < _to_append.getWidth(); i++) {
            for (int j = 0; j < _to_append.getHeight(); j++) {
                if (base.get(i, j)>0 || _to_append.get(i, j)>0){
                    base.set(i, j, 255);
                }
            }
        }
    }

    private static void appendDirections(
            ArrayList[][] direction_map,
            ByteProcessor _to_append,
            int[][] _xy2i,
            int[][] _i2xy,
            int[][] _peaks_i,
            float[][] _branch_score,
            float output_membership_th)
    {

        for (int x=0; x<_to_append.getWidth(); x++) {  // loop elements of the region again
            for (int y = 0; y < _to_append.getHeight(); y++) {

                if (_to_append.get(x, y)>0) {

                    int icoord = _xy2i[x][y];

                    if (icoord!= -1) {  // _peaks_i[icoord] is not null, delineation exists

                        // element of the region is in foreground, take its thetas (if they exist)
                        for (int peak_idx = 0; peak_idx < _peaks_i[icoord].length; peak_idx++) {

                            // check if it exists and if it exists check whether the branch is on
                            int curr_peak_i = _peaks_i[icoord][peak_idx];
                            boolean curr_peak_on = _branch_score[icoord][peak_idx]>output_membership_th; // true;

                            if (curr_peak_i!=-1 && curr_peak_i!=-2 && curr_peak_on) { // indexes of the spatial locations corresponding to peaks

                                int peak_x = _i2xy[curr_peak_i][0]; // PeakExtractor stores spatial location of the follow-up points
                                int peak_y = _i2xy[curr_peak_i][1];

                                if (direction_map[x][y] == null) {
                                    direction_map[x][y] = new ArrayList<float[]>(20);
                                    direction_map[x][y].add(new float[]{peak_x, peak_y});
                                }
                                else {
                                    direction_map[x][y].add(new float[]{peak_x, peak_y});
                                }

                            }

                        }

                    }

                }
            }

        }

    }

    public Overlay[] saveDetection() // will save the detection: ArrayList<CritpointRegion> -> file.det
    {
        // saves detections in predefined output folders
        // 2 outputs:
        // 1. image_name.det textual file with the description of critpoint regions
        //      format: x, y, radius, score, type{BIF,END,CROSS}, dir{vx,vy; vx,vy; ...}
        // 2. return Overlay array that was extracted and saved for each dsens[i] threshold

        Overlay[] ovs = new Overlay[dsens.length]; // alloc output

        for (int i = 0; i < dsens.length; i++) { // loop list of scores with different thresholds

            //// initialize file
            String det_path = output_dir_name.get(i) + File.separator + image_name+".det";

            PrintWriter     logWriter = null;
            try {
                logWriter = new PrintWriter(det_path);
                logWriter.print("");
                logWriter.close();
            } catch (FileNotFoundException ex) {}

            try {
                logWriter = new PrintWriter(new BufferedWriter(new FileWriter(det_path, true)));
                logWriter.println("# DETECTOR...\n# X, Y, RADIUS, SCORE, TYPE, DIRECTIONS(vx,vy, vx,vy...) ");
            } catch (IOException e) {}

            //// initialize Overlay component
            ovs[i] = new Overlay();

            // write the regions
            for (int ii=0; ii<detected_regions.get(i).size(); ii++) {

                if (detected_regions.get(i).get(ii)!=null) { // just in case if the CritpointRegion was null

                    // write one component of the Overlay
                    float cx = detected_regions.get(i).get(ii).centroid[0];
                    float cy = detected_regions.get(i).get(ii).centroid[1];
                    float cr = D[0]/2f;//detected_regions.get(i).get(ii).radius;
                    float sc = detected_regions.get(i).get(ii).score;
                    CritpointRegion.RegionType ctype = detected_regions.get(i).get(ii).type;

                    Color region_color = null;
                    sc = 0.5f;
                    switch (ctype) {
                        case BIF:
                            region_color = new Color(1, 0, 0, sc); // Color.RED;//
                            break;

                        case END:
                            region_color = new Color(1, 1, 0, sc); // Color.YELLOW;//
                            break;

                        case CROSS:
                            region_color = new Color(1, 0, 0, sc); // Color.GREEN;// turned to RED
                            break;

                        default:
                            IJ.log("non valid critical point");
                            break;
                    }

                    if (ctype==null) continue; // skip adding the overlay and line in .det

                    // write line to .det for each region
                    String curr_detection = "";

                    //add region as OvalRoi to output array of Overlays[]
                    OvalRoi ovroi = new OvalRoi(cx-cr+.5, cy-cr+.5, 2*cr, 2*cr);
                    ovroi.setStrokeWidth(1);
                    ovroi.setStrokeColor(region_color);
                    if (true) ovroi.setFillColor(region_color); // ctype== CritpointRegion.RegionType.END
                    ovs[i].add(ovroi);

                    // add line to .det  (what's loaded so far) // HERE IS HOW IT IS WRITTTEN!!!
                    curr_detection +=
                            IJ.d2s(cx,2)+", "+IJ.d2s(cy,2)+", "+
                                    IJ.d2s(cr,2)+", "+
                                    IJ.d2s(sc,2)+", "+detected_regions.get(i).get(ii).type+", ";

                    // add directions (outward_directions) os Line  Overlay component
                    float scale_direction = 2f; // just for visualization
                    int nr_directions = detected_regions.get(i).get(ii).outward_directions.length;
                    for (int j = 0; j < nr_directions; j++) {

                        float dx = detected_regions.get(i).get(ii).outward_directions[j][0];
                        float dy = detected_regions.get(i).get(ii).outward_directions[j][1];

                        curr_detection += IJ.d2s(dx,2)+", "+IJ.d2s(dy,2); // line for .det file
                        if (j<nr_directions-1) curr_detection += ", ";

                        if (!Float.isNaN(dx) && !Float.isNaN(dy)) {

                            Line l = new Line(cx+cr*dx+.5, cy+cr*dy+.5, cx+scale_direction*cr*dx+.5, cy+scale_direction*cr*dy+.5);
                            l.setStrokeWidth(1);
                            l.setStrokeColor(region_color);
                            l.setFillColor(region_color);
                            ovs[i].add(l);

                        }

                    }

                    logWriter.println(curr_detection);

                }

            }

            logWriter.close();
            IJ.log("export " + det_path);

        }

        return ovs;

    }

    public static ReadDET loadDetection(String det_file_path)
    {
        //
        ReadDET read_det = new ReadDET(det_file_path);

//        ArrayList<CritpointRegion> regs = new ArrayList<CritpointRegion>();

        return read_det;

    }

    private void appendDetectedRegions( // appends to the list of CritpointRegion, without directions
                                        ByteProcessor               _region_map,
                                        FloatProcessor              _score_map,
                                        CritpointRegion.RegionType  _choose_type,
                                        ArrayList<CritpointRegion>  region_list
    ) {

        // check the number of pixels that were above the threshold
        byte[] check_vals = (byte[])_region_map.getPixels();
        int cnt = 0;
        for (int i = 0; i < check_vals.length; i++) if ((check_vals[i]&0xff) > 0) cnt++;
        if (cnt/(float)check_vals.length>ratioOnPix) {
            IJ.log("break! overshot, %ONpix = " + (cnt/(float)check_vals.length));
            return;
        }

        ///// connected components
        ip_exporter.setProcessor(_region_map);
        Conn conn_reg = new Conn(ip_exporter, true);
        conn_reg.run("");
        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

        int Dmax = (int) Math.round(Stat.get_max(D));
        int Amax = (int) (2*Math.pow(Dmax, 2)); // 2 is heuristic
//		int Amin = 4;

//		if (regs.size()>maxNrRegions) {
//			IJ.log("break! overshot, #regions = " + regs.size());
//			return;
//		}

        // regs
        for (int i=0; i<regs.size(); i++) {

//            if (regs.get(i).size() > Amax) continue; // go to the next region
            if (regs.get(i).size() < 2) continue; // go to the next region if it's only 2 pixels

            float Cx=0, Cy=0; 	// centroid
            float C=0;        	// score
            float Cr; 			// radius

            CritpointRegion.RegionType Ctype;   	// type
            float[][] Cdirections;           		// directions

            for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region

                int xcoord = regs.get(i).get(aa)[1];
                int ycoord = regs.get(i).get(aa)[0];

                Cx += regs.get(i).get(aa)[1];
                Cy += regs.get(i).get(aa)[0];
                C += _score_map.getf(xcoord, ycoord); //_critpoint_det[xcoord][ycoord];

            }

            Cx /= regs.get(i).size();
            Cy /= regs.get(i).size();   // centroid

            C /= regs.get(i).size();    // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)

            Cr = (float) Math.ceil(             1f*Math.sqrt(regs.get(i).size()/3.14f)        ); // radius is wrt to the area
            Cr = (Cr<1)? 1 : Cr;

            if (_choose_type == CritpointRegion.RegionType.END) {

                Ctype = CritpointRegion.RegionType.END;
                Cdirections = new float[1][2];
                Arrays.fill(Cdirections[0], Float.NaN);

                region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 1)); // take one from Cdirections

            }
            else if (_choose_type == CritpointRegion.RegionType.BIF_CROSS) {

                Ctype = CritpointRegion.RegionType.BIF;
                Cdirections = new float[3][2];
                Arrays.fill(Cdirections[0], Float.NaN);
                Arrays.fill(Cdirections[1], Float.NaN);
                Arrays.fill(Cdirections[2], Float.NaN);

                region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 3)); // take three from Cdirections

            }
        }
    }

    private void appendDetectedRegions( // appends to the list of CritpointRegion, extract directions
                                        ByteProcessor               _region_map,
                                        ArrayList[][]               _peaks_map,
                                        FloatProcessor              _score_map,
                                        CritpointRegion.RegionType  _choose_type,
                                        float                       alfa_deg,       // mean-shift param (for directionality analysis)
                                        ArrayList<CritpointRegion>  region_list     // destination to append
    )
    {

        ///// connected components
        ip_exporter.setProcessor(_region_map);
        Conn conn_reg = new Conn(ip_exporter, true);
        conn_reg.run("");
        ArrayList<ArrayList<int[]>> regs = conn_reg.getConnectedRegions();

        int Dmax = (int) Math.round(Stat.get_max(D));
        int Amax = (int) (2*Math.pow(Dmax, 2)); // 2 is heuristic
        int nr_regs = 0;

        // regs
        for (int i=0; i<regs.size(); i++) {

            if (regs.get(i).size() > Amax) continue; // go to the next region
            if (regs.get(i).size() < 4) continue; // go to the next region

            nr_regs++;

            float Cx=0, Cy=0; 	// centroid
            float C=0;        	// score
            float Cr; 			// radius

            CritpointRegion.RegionType Ctype;   	// type
            float[][] Cdirections;           		// directions

            for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region

                int xcoord = regs.get(i).get(aa)[1];
                int ycoord = regs.get(i).get(aa)[0];

                Cx += regs.get(i).get(aa)[1];
                Cy += regs.get(i).get(aa)[0];
                C += _score_map.getf(xcoord, ycoord); //_critpoint_det[xcoord][ycoord];

            }

            Cx /= regs.get(i).size();
            Cy /= regs.get(i).size();   // centroid

            C /= regs.get(i).size();    // score (calculate it here regardless of the type, says how much average fuzzy score was confident on the output)

            Cr = (float) Math.ceil(  (1f * Math.sqrt(regs.get(i).size()/3.14f))  ); // radius is wrt to the area
            Cr = (Cr<1)? 1 : Cr;

            Ctype = null;   			// values to add for this region

            ArrayList<float[]> vxy = new ArrayList<float[]>();
            vxy.clear();    // list of local directions taken from the region, clear before starting for each region

            for (int aa=0; aa<regs.get(i).size(); aa++) {  // loop elements of the region again to combine the directions together

                // every location will have list of peaks (2d locations surrounding the point)
                int xcoord = regs.get(i).get(aa)[1];
                int ycoord = regs.get(i).get(aa)[0];

                if (_peaks_map[xcoord][ycoord]!=null) {

                    for (int j = 0; j < _peaks_map[xcoord][ycoord].size(); j++) {

                        float[] take_xy = (float[]) _peaks_map[xcoord][ycoord].get(j);    // read it from the map

                        float peak_x = take_xy[0];
                        float peak_y = take_xy[1];

                        float[] unit_vxy = new float[]{peak_x-Cx, peak_y-Cy};
                        float norm_vxy = (float) Math.sqrt(Math.pow(unit_vxy[0],2)+Math.pow(unit_vxy[1],2));
                        unit_vxy[0] = unit_vxy[0] / norm_vxy;
                        unit_vxy[1] = unit_vxy[1] / norm_vxy;

                        vxy.add(unit_vxy); // vxy.add(new float[]{take_xy[0], take_xy[1]});

                    }

                }

            }

            // vxy list is formed.. to make sense - take those that had more than certain amount of directions
            // to be sure that the mean shift makes sense and the region itself is salient to be added at all
            if (_choose_type == CritpointRegion.RegionType.END && vxy.size()>0) {

                ArrayList<float[]> direction_clusters_vxy_count = ClusterDirections.run(vxy, alfa_deg);
                Ctype = CritpointRegion.RegionType.END; // this is known

                Cdirections = new float[direction_clusters_vxy_count.size()][2];  // will take 2 columns (3rd will be used to determine if it is crossing)
                for (int k=0; k<direction_clusters_vxy_count.size(); k++) {
                    Cdirections[k][0] = direction_clusters_vxy_count.get(k)[0];
                    Cdirections[k][1] = direction_clusters_vxy_count.get(k)[1];
                }

                region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 1)); // take one from the top only, no need to discriminate further

            }
            else if (_choose_type == CritpointRegion.RegionType.BIF_CROSS && vxy.size()>2) {

                ArrayList<float[]> direction_clusters_vxy_count = ClusterDirections.run(vxy, alfa_deg);

                //direction_clusters_vxy_count.get( RANKING )[X_Y_COUNT]
                // decide whether it is bif (3- clusters) or CROSS (4 clusters)
                if (direction_clusters_vxy_count.size()==4) { // if the last one is balanced towards smallest remaining

                    boolean is_cross = direction_clusters_vxy_count.get(3)[2] / direction_clusters_vxy_count.get(2)[2] >= 0.8;

                    if (is_cross) {
                        Ctype = CritpointRegion.RegionType.CROSS; // (add all 4)
                    } else {
                        Ctype = CritpointRegion.RegionType.BIF;   // (add all 3)
                    }

                }
                else if (direction_clusters_vxy_count.size()==3) {
                    Ctype = CritpointRegion.RegionType.BIF;
                }
                else {
                    // region was classified as junction (bif. or crs.) but the delineation directions were not 3 at least after mean-shifting
                    Ctype = CritpointRegion.RegionType.BDY;   // it's not END, it's not BIF, or CROSS
                }

                Cdirections = new float[direction_clusters_vxy_count.size()][2];  // will take 2 columns (3rd will be used to determine if it is crossing)
                for (int k=0; k<direction_clusters_vxy_count.size(); k++) {
                    Cdirections[k][0] = direction_clusters_vxy_count.get(k)[0];
                    Cdirections[k][1] = direction_clusters_vxy_count.get(k)[1];
                }

                switch (Ctype) {
                    case BIF:
                        region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 3));
                        break;
                    case CROSS:
                        region_list.add(new CritpointRegion(Ctype, Cx, Cy, Cr, C, Cdirections, 4));
                        break;
                    case BDY:
                        // don't add it
                        break;
                    default:
                        break;
                }

            } // else nothing

        }

        IJ.log("[" + nr_regs + "/" + regs.size() + " conn. regions taken]");

    }

    public static void print(int atX, int atY)
    {

        int atLoc = Masker.xy2i[atX][atY];

        IJ.log(String.format("/**** DETECTOR 2D TOOL (%5d, %5d) [%10d] ****/", atX, atY, atLoc));

        if (atLoc != -1) {  // first level check, intensity masker

            if (Delineator.delin2[atLoc]!=null) {  // second level check, profile masker

                String printout = "";

                printout += "\nDELINEATION INDEXES:\n";

                for (int ii = 0; ii< Delineator.delin2[atLoc].length; ii++) {
                    printout += ii+"\t->\t";
                    for (int jj = 0; jj< Delineator.delin2[atLoc][ii].length; jj++) {

                        if (Delineator.delin2[atLoc][ii][jj]==-1) {
                            printout += "BGRD";
                        }
                        else if (Delineator.delin2[atLoc][ii][jj]==-2) {
                            printout += "NONE";
                        }
                        else {
                            printout += IJ.d2s(Delineator.delin2[atLoc][ii][jj], 2);
                        }

                        if (jj== Delineator.delin2[atLoc][ii].length-1) printout += "\n";
                        else printout += ",  ";
                    }
                }

                printout += "\nREFINED LOCS:\n";
                if (Delineator.xy2[atLoc]!=null) {

                    for (int b = 0; b< Delineator.xy2[atLoc].length; b++) {

                        printout += b+"\t->\t";

                        if (Delineator.xy2[atLoc][b]!=null) {

                            if (Delineator.xy2[atLoc][b][0] != null) {

                                for (int l = 0; l< Delineator.xy2[atLoc][b][0].length; l++) {

                                    printout += "("+IJ.d2s(Delineator.xy2[atLoc][b][0][l], 2)+", "+IJ.d2s(Delineator.xy2[atLoc][b][1][l], 2)+")";

                                    if (l== Delineator.xy2[atLoc][b][0].length-1) printout += "\n";
                                    else printout += ", ";
                                }

                            }
                            else {

                                printout += "NULL\n";

                            }




                        }
                        else {
                            printout += "NONE\n";
                        }

                    }

                }
                else {
                    printout += "SKIPPED CALCULATING HERE (THERE WAS A THREAD POINTING TO BGRD)\n";
                }

                printout += "\nREFINED VECS:\n";
                if (Delineator.vxy2[atLoc]!=null) {

                    for (int b = 0; b< Delineator.vxy2[atLoc].length; b++) {

                        printout += b+"\t->\t";

                        if (Delineator.vxy2[atLoc][b]!=null) {

                            if (Delineator.vxy2[atLoc][b][0] != null) {

                                for (int l = 0; l< Delineator.vxy2[atLoc][b][0].length; l++) {

                                    printout += "("+IJ.d2s(Delineator.vxy2[atLoc][b][0][l], 2)+", "+IJ.d2s(Delineator.vxy2[atLoc][b][1][l], 2)+")";

                                    if (l== Delineator.vxy2[atLoc][b][0].length-1) printout += "\n";
                                    else printout += ", ";
                                }

                            }
                            else {

                                printout += "NULL\n";

                            }


                        }
                        else {
                            printout += "NONE\n";
                        }

                    }

                }
                else {
                    printout += "SKIPPED CALCULATING HERE (THERE WAS A THREAD POINTING TO BGRD)\n";
                }

                printout += "\nLHOOD STHNESS NCC\n";

                if (Delineator.smoothness[atLoc]!=null && Ncc.scores[atLoc]!=null) {
                    for (int b=0; b<4; b++)
                        printout += b + " -> " + IJ.d2s(PeakExtractor.peaks_lhood[atLoc][b], 1) + " " + IJ.d2s(Delineator.smoothness[atLoc][b], 1) + " " + IJ.d2s(Ncc.scores[atLoc][b], 1) + "\n";
                }
                else {
                    printout += "NONE(profile was flat)\n";
                }

                printout += "\nEND <- NONE -> JUN\n";

                printout += IJ.d2s(FuzzyDetector.endpoint_score[atLoc],2) + "<- NONE ->" + IJ.d2s(FuzzyDetector.junction_score[atLoc],2);

//                if (FuzzyDetector.endpoint_score[atLoc]!=null && Ncc.scores[atLoc]!=null) {
//                    for (int b=0; b<4; b++)
//                        printout += b + " -> " + IJ.d2s(PeakExtractor.peaks_lhood[atLoc][b], 1) + " " + IJ.d2s(Delineator.smoothness[atLoc][b], 1) + " " + IJ.d2s(Ncc.scores[atLoc][b], 1) + "\n";
//                }
//                else {
//                    printout += "NONE(profile was flat)\n";
//                }

                IJ.log(printout);

            }
            else {
                IJ.log("not enough variation in profile, profile masker said so");
            }

        }
        else {
            IJ.log("background point, intensity masker said so");
        }

    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1)
    {
        // emulate what the Fuzzy instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(ncc_1, likelihood_1, smoothness_1, dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2)
    {
        // emulate what the Fuzzy instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2,
                                             float ncc_3, float likelihood_3, float smoothness_3)
    {
        // emulate what the Fuzzy instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];
        float[] dummy1 = new float[4];
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                ncc_3, likelihood_3, smoothness_3,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    private ImageStack getFuzzyLogicResponse(float ncc_1, float likelihood_1, float smoothness_1,
                                             float ncc_2, float likelihood_2, float smoothness_2,
                                             float ncc_3, float likelihood_3, float smoothness_3,
                                             float ncc_4, float likelihood_4, float smoothness_4)
    {
        // emulate what the Fuzzy instance was giving for this input, use sample_fls initiated in the same way
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = true; // enable capturing the results of different stages of fuzzy logic detection
        float[] dummy = new float[3];           // end, none, jun score
        float[] dummy1 = new float[4];          //
        sample_fls.critpointScore(
                ncc_1, likelihood_1, smoothness_1,
                ncc_2, likelihood_2, smoothness_2,
                ncc_3, likelihood_3, smoothness_3,
                ncc_4, likelihood_4, smoothness_4,
                dummy, dummy1);
        ImageStack is_out; // retrieve the cummulative log
        is_out = sample_fls.fls_steps.duplicate();
        sample_fls.clearLog(); // to be empty for the next call
        sample_fls.verbose = false;
        return is_out;
    }

    public ImageStack getFuzzyLogicResponse(int atX, int atY)
    {

        ImageStack is_out = new ImageStack(plotw, ploth);

        int atLoc = Masker.xy2i[atX][atY];

        if (atLoc != -1) {

            if (Delineator.smoothness[atLoc]!=null && Ncc.scores[atLoc]!=null) {

                int cnt = 0; // count branches that are existing
                float ncc_1=0, lhood_1=0, smthness_1=0;
                float ncc_2=0, lhood_2=0, smthness_2=0;
                float ncc_3=0, lhood_3=0, smthness_3=0;
                float ncc_4=0, lhood_4=0, smthness_4=0;

                for (int b=0; b<4; b++) {

                    float curr_ncc      = Ncc.scores[atLoc][b];
                    float curr_lhood    = PeakExtractor.peaks_lhood[atLoc][b];
                    float curr_smooth   = Delineator.smoothness[atLoc][b];

                    if (!Float.isNaN(curr_ncc) && !Float.isNaN(curr_lhood) && !Float.isNaN(curr_smooth)) {

                        cnt++;

                        if (cnt==1) {
                            ncc_1 = curr_ncc;
                            lhood_1 = curr_lhood;
                            smthness_1 = curr_smooth;
                        }
                        else if (cnt==2) {
                            ncc_2 = curr_ncc;
                            lhood_2 = curr_lhood;
                            smthness_2 = curr_smooth;
                        }
                        else if (cnt==3) {
                            ncc_3 = curr_ncc;
                            lhood_3 = curr_lhood;
                            smthness_3 = curr_smooth;
                        }
                        else if (cnt==4) {
                            ncc_4 = curr_ncc;
                            lhood_4 = curr_lhood;
                            smthness_4 = curr_smooth;
                        }


                    }

                }

                if (cnt==1) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1);
                }
                else if (cnt==2) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2);
                }
                else if (cnt==3) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3);
                }
                else if (cnt==4) {
                    is_out = getFuzzyLogicResponse(ncc_1, lhood_1, smthness_1, ncc_2, lhood_2, smthness_2, ncc_3, lhood_3, smthness_3, ncc_4, lhood_4, smthness_4);
                }

            }
            else {
                is_out.addSlice("NOTHING", new Plot("","","").getProcessor());
            }

        }
        else {
            is_out.addSlice("NOTHING", new Plot("","","").getProcessor());
        }

        return is_out;

    }

}
