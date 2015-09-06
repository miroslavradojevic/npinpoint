import ij.*;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

/**
 * Created by miroslav on 5-9-15.
 */
public class Critpoint2D implements PlugIn, MouseListener, MouseMotionListener {

    String      image_path;

    /*
    interface elements - all the windows that pop up as you click/move with mouse
     */
    ImagePlus       pfl_im  = new ImagePlus();
    ImagePlus       pfl_im1 = new ImagePlus();
    ImagePlus       pfl_im2 = new ImagePlus();
    ImagePlus       pfl_im3 = new ImagePlus();
    ImagePlus       pfl_im4 = new ImagePlus();

    ImageStack      pfl_is  = null;
    ImageStack      pfl_is1 = null;
    ImageStack      pfl_is2 = null;
    ImageStack      pfl_is3 = null;
    ImageStack      pfl_is4 = null;

    ImageWindow     pfl_iw4 = null;
    ImageWindow   	pfl_iw3 = null;
    ImageWindow   	pfl_iw2 = null;

    boolean first4 = true;
    boolean first3 = true;
    boolean first2 = true;

    ImageCanvas cnv;

    int plot_w = 528;
    int plot_h = 255;
    int upper_left_x = 70;
    int upper_left_y = 50;

    Detector2D det2d = null;

    public void run(String sss){

        // load the image through the menu
        String in_folder = Prefs.get("id.folder", System.getProperty("user.home"));
        OpenDialog.setDefaultDirectory(in_folder);
        OpenDialog dc = new OpenDialog("Select file");
        in_folder = dc.getDirectory();
        image_path = dc.getPath();
        if (image_path==null) return;
        Prefs.set("id.folder", in_folder);

        ImagePlus ip_load = new ImagePlus(image_path);
        if(ip_load==null) return;

        /*
        	load the detection parameters
         */
        boolean _show_junctions, _show_endpoints, _enable_interactive, _save_midresults, _calculate_directions;
        float _s;
        float _sigma_ratio;
        String _Dlist       = "";
        String _dsens_list  = "";
        float _ncc_high, _ncc_low, _likelihood_high, _likelihood_low, _smoothness_low, _smoothness_high;
//        float _detection_sensitivity; // ranges from 0 (max) - 1 (max. entropy) will be used at thresholding

        // check if the parameters were submitted through the macro before rising up the Generic Dialog
        // enables calling plugin from the macro without opening the graphical window
        // useful to call plugins with parameters submitted by ij macro in fiji headless mode
        if (Macro.getOptions()==null) {

            // generic dialog (graphic)
            _show_junctions 		= 			Prefs.get("critpoint.detection2d.show_junctions", true);
            _show_endpoints 		= 			Prefs.get("critpoint.detection2d.show_endpoints", true);
            _enable_interactive     = 			Prefs.get("critpoint.detection2d.enable_interactive", true);
            _calculate_directions   =           Prefs.get("critpoint.detection2d.calculate_directions", false);
            _save_midresults        =           Prefs.get("critpoint.detection2d.save_midresults", false);
            _s						= (float)	Prefs.get("critpoint.detection2d.s", 1.2f);
            _Dlist 					= 			Prefs.get("critpoint.detection2d.d", "6");
            _sigma_ratio			= (float) 	Prefs.get("critpoint.detection2d.sigma_ratio", 		1f/6f);
            _ncc_high               = (float)   Prefs.get("critpoint.detection2d.ncc_high", 		0.90f);
            _ncc_low				= (float) 	Prefs.get("critpoint.detection2d.ncc_low",	 		0.50f);
            _likelihood_high        = (float)   Prefs.get("critpoint.detection2d.likelihood_high", 	0.90f);
            _likelihood_low			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low", 	0.30f);
            _smoothness_high        = (float)   Prefs.get("critpoint.detection2d.smoothness_high", 	20f);
            _smoothness_low         = (float)   Prefs.get("critpoint.detection2d.smoothness_low", 	0f);
            _dsens_list             =           Prefs.get("critpoint.detection2d.detection_sensitivity", "1.0"); // start with the highest sensitivity

            GenericDialog gd = new GenericDialog("DETECTOR2D");

            gd.addCheckbox("JUNCTIONS", 		    _show_junctions);
            gd.addCheckbox("ENDPOINTS", 		    _show_endpoints);
            gd.addCheckbox("INTERACTIVE",           _enable_interactive);
            gd.addCheckbox("CALCULATE_DIRECTIONS",  _calculate_directions);
            gd.addCheckbox("SAVE_MIDRESULTS",       _save_midresults);

            gd.addStringField("Dlist", 				_Dlist, 30);
            gd.addNumericField("sigma_ratio",       _sigma_ratio, 	    2,  10, "");
            gd.addNumericField("s", 				_s,					1,	10,	"");
            gd.addNumericField("NCC_HIGH", 	        _ncc_high, 			2,  10, "");
            gd.addNumericField("NCC_LOW",           _ncc_low, 		    2,  10, "");
            gd.addNumericField("LIKELIHOOD_HIGH", 	_likelihood_high, 	2,  10, "");
            gd.addNumericField("LIKELIHOOD_LOW",    _likelihood_low, 	2,  10, "");
            gd.addNumericField("SMOOTHNESS_HIGH", 	_smoothness_high, 	2,  10, "");
            gd.addNumericField("SMOOTHNESS_LOW",    _smoothness_low, 	2,  10, "");
            gd.addStringField("dsens_list", _dsens_list, 30);//addNumericField("DETECTION SENSITIVITY", _detection_sensitivity, 2, 10, "");

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _show_junctions = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_junctions", 		_show_junctions);
            _show_endpoints = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_endpoints", 		_show_endpoints);
            _enable_interactive = gd.getNextBoolean();				Prefs.set("critpoint.detection2d.enable_interactive", 	_enable_interactive);
            _calculate_directions =  gd.getNextBoolean();			Prefs.set("critpoint.detection2d.calculate_directions", _calculate_directions);
            _save_midresults = gd.getNextBoolean();					Prefs.set("critpoint.detection2d.save_midresults", 		_save_midresults);

            _Dlist       	= gd.getNextString(); 				    Prefs.set("critpoint.detection2d.d", _Dlist);
            _sigma_ratio	= (float) gd.getNextNumber();			Prefs.set("critpoint.detection2d.sigma_ratio", _sigma_ratio);
            _s          	= (float) gd.getNextNumber(); 			Prefs.set("critpoint.detection2d.s", _s);
            _ncc_high   	= (float) gd.getNextNumber();   	    Prefs.set("critpoint.detection2d.ncc_high", _ncc_high);
            _ncc_low    	= (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low", _ncc_low);
            _likelihood_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_high", 	_likelihood_high);
            _likelihood_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_low", 	_likelihood_low);
            _smoothness_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_high", 	_smoothness_high);
            _smoothness_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_low", 	_smoothness_low);
            _dsens_list     = gd.getNextString();	                Prefs.set("critpoint.detection2d.detection_sensitivity", _dsens_list);

        }
        else { // continue with macro arguments without rising graphic window

//            System.out.println(Macro.getOptions());

            _show_junctions = Boolean.valueOf(Macro.getValue(Macro.getOptions(),        "junctions", String.valueOf(false)));
            _show_endpoints = Boolean.valueOf(Macro.getValue(Macro.getOptions(),        "endpoints", String.valueOf(false)));
            _enable_interactive = Boolean.valueOf(Macro.getValue(Macro.getOptions(),    "interactive", String.valueOf(false)));
            _calculate_directions = Boolean.valueOf(Macro.getValue(Macro.getOptions(),  "calculate_directions", String.valueOf(false)));
            _save_midresults = Boolean.valueOf(Macro.getValue(Macro.getOptions(),       "save_midresults", String.valueOf(false)));
            _Dlist = Macro.getValue(Macro.getOptions(),                                 "dlist", String.valueOf(6));
            _sigma_ratio = Float.valueOf(Macro.getValue(Macro.getOptions(),             "sigma_ratio", String.valueOf(0.25)));
            _s = Float.valueOf(Macro.getValue(Macro.getOptions(),                       "s", String.valueOf(1.1)));

            _ncc_high           = Float.valueOf(Macro.getValue(Macro.getOptions(),      "ncc_high", String.valueOf(0.95)));
            _ncc_low            = Float.valueOf(Macro.getValue(Macro.getOptions(),      "ncc_low", String.valueOf(0.2)));
            _likelihood_high    = Float.valueOf(Macro.getValue(Macro.getOptions(),      "likelihood_high", String.valueOf(0.5)));
            _likelihood_low     = Float.valueOf(Macro.getValue(Macro.getOptions(),      "likelihood_low", String.valueOf(0.0)));
            _smoothness_high    = Float.valueOf(Macro.getValue(Macro.getOptions(),      "smoothness_high", String.valueOf(10)));
            _smoothness_low     = Float.valueOf(Macro.getValue(Macro.getOptions(),      "smoothness_low", String.valueOf(20)));
            _dsens_list         = Macro.getValue(Macro.getOptions(),                    "dsens_list", String.valueOf(1.0)); // only one value by default

        }

//        System.out.println("----------");
//        System.out.println("loaded:");
//        System.out.println(_dsens_list);
//        System.out.println(_Dlist);
//        System.out.println("----------");

        String[] dd = _Dlist.split(","); if (dd.length==0) return;
        float[] _D = new float[dd.length];
        for (int i=0; i<dd.length; i++) _D[i] = Float.valueOf(dd[i]);

        dd = _dsens_list.split(",");  if (dd.length==0) return;
        float[] _dsens = new float[dd.length];
        for (int i = 0; i < dd.length; i++) _dsens[i] = Float.valueOf(dd[i]);

        // detection
        System.out.print("\nDetector2D... ");
        det2d = new Detector2D(
                ip_load,
                _s,
                _D,
                _sigma_ratio,
                _ncc_high,
                _ncc_low,
                _likelihood_high,
                _likelihood_low,
                _smoothness_high,
                _smoothness_low ,
                _dsens
        );

        // they are initialized by default already , this just overwrites
        det2d.saveMidresults(_save_midresults); 	 // t
        det2d.do_directions = _calculate_directions;    // f
        det2d.do_endpoints = _show_endpoints;        // t
        det2d.do_junctions = _show_junctions;        // t

        det2d.run();
        Overlay[] detection_overlays = det2d.saveDetection();  // export .zip and .det with detections

        if (_enable_interactive) {

            // allocate image stack
            ImageStack out_stack = new ImageStack(ip_load.getWidth(), ip_load.getHeight());
            // allocate joint Overlay for the stack
            Overlay    out_ovl = new Overlay();

            for (int i = detection_overlays.length-1; i >= 0; i--) {  // loops dsens.length times

                if (detection_overlays[i].size()>0) {

                    out_stack.addSlice("DET2D,dsens="+IJ.d2s(_dsens[i],2)+".tif", ip_load.getProcessor().duplicate());

                    for (int j = 0; j < detection_overlays[i].size(); j++) {

                        Roi get_roi = detection_overlays[i].get(j);
                        get_roi.setPosition(out_stack.getSize());
                        out_ovl.add(get_roi);

                    }

                }
                else System.out.print("Empty overlay with detections.");

            }

            ImagePlus out_imp = new ImagePlus("Critpoint2D", out_stack);
            out_imp.setOverlay(out_ovl);
            out_imp.show();
//            out_imp.getCanvas().zoomIn(0,0);
//            out_imp.getCanvas().zoomIn(0,0);
//            out_imp.getCanvas().zoomIn(0,0);
//            out_imp.getCanvas().zoomIn(0,0);
//            out_imp.getWindow().maximize();

            ip_load.show();
            cnv = ip_load.getCanvas();
            cnv.addMouseListener(this);
            cnv.addMouseMotionListener(this);
            IJ.setTool("hand");
        }
        else { // if not interactive just show the stack with detections

            // allocate image stack

            ImageStack out_stack = new ImageStack(ip_load.getWidth(), ip_load.getHeight());

            Overlay    out_ovl = new Overlay();

            for (int i = detection_overlays.length-1; i >= 0; i--) {  // loops dsens.length times

                if (detection_overlays[i].size()>0) {

                    out_stack.addSlice("DET2D,dsens="+ IJ.d2s(_dsens[i], 2)+".tif", ip_load.getProcessor().duplicate());

                    for (int j = 0; j < detection_overlays[i].size(); j++) {

                        Roi get_roi = detection_overlays[i].get(j);
                        get_roi.setPosition(out_stack.getSize());
                        out_ovl.add(get_roi);

                    }

                }
                else System.out.print("Empty overlay with detections.");

            }

            ImagePlus out_imp = new ImagePlus("Critpoint2D", out_stack);
            out_imp.setOverlay(out_ovl);
            out_imp.show();

        }

        System.out.println("DONE.");

    }

    public void mouseClicked(MouseEvent e)
    {
//        mouseMoved(e);
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());


        pfl_is3 = Delineator2D.getClusterMask(clickX, clickY);
        if (pfl_is3.getSize()>0) {
            pfl_im3.setStack("cluster_mask", pfl_is3);
            if (first3) {
                pfl_im3.show();
                pfl_iw3 = pfl_im3.getWindow();
                pfl_iw3.setLocation(upper_left_x+plot_w+20, upper_left_y);
                first3 = false;
            }
            pfl_im3.setStack("fitting_scores", pfl_is3);
            //ImageWindow iw = pfl_im3.getWindow();
            pfl_im3.updateAndDraw();
//            pfl_im3.show();
        }

        pfl_is2 = (det2d!=null)? det2d.getFuzzyLogicResponse(clickX, clickY) : null;
        if (pfl_is2!=null && pfl_is2.getSize()>0) {
            pfl_im2.setStack("what happened in fuzzy", pfl_is2);
            if (first2) {
                pfl_im2.show();
                pfl_iw2 = pfl_im2.getWindow();
                pfl_iw2.setLocation(upper_left_x, upper_left_y+plot_h+50);
                first2 = false;
            }
            pfl_im2.setStack("what happened in fuzzy", pfl_is2);
//            ImageWindow iw = pfl_im2.getWindow();
            pfl_im2.updateAndDraw();
        }

        IJ.setTool("hand");

    }

    public void mouseMoved(MouseEvent e)
    {
        //mouseClicked(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        /*
            output Overlay & update canvas with the original
         */
        Overlay ov_to_add = Delineator2D.getDelineationOverlay(clickX, clickY);
        cnv.setOverlay(ov_to_add);

        Detector2D.print(clickX, clickY); // print extracted features

        pfl_is4 = PeakExtractor2D.getProfileWithPeaks(clickX, clickY);
        if (pfl_is4.getSize()>0) pfl_im4.setStack("local_profile_with_peaks", pfl_is4);
        if (first4) {
            pfl_im4.show();
            pfl_iw4 = pfl_im4.getWindow();
            pfl_iw4.setLocation(upper_left_x, upper_left_y);
            first4 = false;
        }

        IJ.setTool("hand");

    }

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}
