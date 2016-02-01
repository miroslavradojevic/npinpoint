package com.braincadet.npinpoint;

import ij.*;
import ij.gui.*;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;

import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;

public class Critpoint implements PlugIn, MouseListener, MouseMotionListener {

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

    int plot_w = new Plot("","","",new float[1],new float[1]).getSize().width;
    int plot_h = new Plot("","","",new float[1],new float[1]).getSize().height;
    int upper_left_x = 20;
    int upper_left_y = 20;

    Detector det2d = null;

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
        if (ip_load.getStack().getSize()>1) {IJ.log("Image needs to be 2D."); return;}

        /*
        	load the detection parameters
         */
        boolean _show_junctions = true, _show_endpoints = true, _enable_interactive = true, _save_midresults = true, _calculate_directions = false;
//        float _s;
//        float _sigma_ratio;
        String _Dlist       = "6";
        String _dsens_list  = "0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1";
        float _ncc_high = 0.95f, _ncc_low = 0.4f, _likelihood_high = 0.4f, _likelihood_low = 0.00f, _smoothness_low = 20f, _smoothness_high = 0f;

        // check if the parameters were submitted through the macro before rising up the Generic Dialog
        // enables calling plugin from the macro without opening the graphical window
        // useful to call plugins with parameters submitted by ij macro in fiji headless mode
        if (Macro.getOptions()==null) {

            // generic dialog (graphic)
            _show_junctions 		= 			Prefs.get("critpoint.detection2d.show_junctions",       _show_junctions);
            _show_endpoints 		= 			Prefs.get("critpoint.detection2d.show_endpoints",       _show_endpoints);
            _enable_interactive     = 			Prefs.get("critpoint.detection2d.enable_interactive",   _enable_interactive);
            _calculate_directions   =           Prefs.get("critpoint.detection2d.calculate_directions", _calculate_directions);
            _save_midresults        =           Prefs.get("critpoint.detection2d.save_midresults",      _save_midresults);
//            _s						= (float)	Prefs.get("critpoint.detection2d.s",                1.2f);
            _Dlist 					= 			Prefs.get("critpoint.detection2d.d",                    _Dlist);
//            _sigma_ratio			= (float) 	Prefs.get("critpoint.detection2d.sigma_ratio", 		    1f/6f);
            _ncc_high               = (float)   Prefs.get("critpoint.detection2d.ncc_high", 		    _ncc_high);
            _ncc_low				= (float) 	Prefs.get("critpoint.detection2d.ncc_low",	 		    _ncc_low);
            _likelihood_high        = (float)   Prefs.get("critpoint.detection2d.likelihood_high", 	    _likelihood_high);
            _likelihood_low			= (float) 	Prefs.get("critpoint.detection2d.likelihood_low", 	    _likelihood_low);
            _smoothness_high        = (float)   Prefs.get("critpoint.detection2d.smoothness_high", 	    _smoothness_high);
            _smoothness_low         = (float)   Prefs.get("critpoint.detection2d.smoothness_low", 	    _smoothness_low);
            _dsens_list             =           Prefs.get("critpoint.detection2d.detection_sensitivity", _dsens_list); // start with the highest sensitivity

            GenericDialog gd = new GenericDialog("DETECTOR");

            gd.addCheckbox("JUNCTIONS", 		    _show_junctions);
            gd.addCheckbox("ENDPOINTS", 		    _show_endpoints);
            gd.addCheckbox("INTERACTIVE",           _enable_interactive);
            gd.addCheckbox("CALCULATE_DIRECTIONS",  _calculate_directions);
            gd.addCheckbox("SAVE_MIDRESULTS",       _save_midresults);

            gd.addStringField("Dlist", 				_Dlist,             20);
//            gd.addNumericField("sigma_ratio",       _sigma_ratio, 	    2,  10, "");
//            gd.addNumericField("s", 				_s,					1,	10,	"");
            gd.addNumericField("NCC_HIGH", 	        _ncc_high, 			2,  10, "");
            gd.addNumericField("NCC_LOW",           _ncc_low, 		    2,  10, "");
            gd.addNumericField("LIKELIHOOD_HIGH", 	_likelihood_high, 	2,  10, "");
            gd.addNumericField("LIKELIHOOD_LOW",    _likelihood_low, 	2,  10, "");
            gd.addNumericField("SMOOTHNESS_HIGH", 	_smoothness_high, 	0,  10, "");
            gd.addNumericField("SMOOTHNESS_LOW",    _smoothness_low, 	0,  10, "");
            gd.addStringField("dsens_list",         _dsens_list,        30);

            gd.showDialog();
            if (gd.wasCanceled()) return;

            _show_junctions = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_junctions", 		_show_junctions);
            _show_endpoints = gd.getNextBoolean();      	    	Prefs.set("critpoint.detection2d.show_endpoints", 		_show_endpoints);
            _enable_interactive = gd.getNextBoolean();				Prefs.set("critpoint.detection2d.enable_interactive", 	_enable_interactive);
            _calculate_directions =  gd.getNextBoolean();			Prefs.set("critpoint.detection2d.calculate_directions", _calculate_directions);
            _save_midresults = gd.getNextBoolean();					Prefs.set("critpoint.detection2d.save_midresults", 		_save_midresults);

            _Dlist       	= gd.getNextString(); 				    Prefs.set("critpoint.detection2d.d", _Dlist);
//            _sigma_ratio	= (float) gd.getNextNumber();			Prefs.set("critpoint.detection2d.sigma_ratio", _sigma_ratio);
//            _s          	= (float) gd.getNextNumber(); 			Prefs.set("critpoint.detection2d.s", _s);
            _ncc_high   	= (float) gd.getNextNumber();   	    Prefs.set("critpoint.detection2d.ncc_high", _ncc_high);
            _ncc_low    	= (float) gd.getNextNumber();  		    Prefs.set("critpoint.detection2d.ncc_low", _ncc_low);
            _likelihood_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_high", 	_likelihood_high);
            _likelihood_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.likelihood_low", 	_likelihood_low);
            _smoothness_high= (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_high", 	_smoothness_high);
            _smoothness_low = (float) gd.getNextNumber();   		Prefs.set("critpoint.detection2d.smoothness_low", 	_smoothness_low);
            _dsens_list     = gd.getNextString();	                Prefs.set("critpoint.detection2d.detection_sensitivity", _dsens_list);

        }
        else { // continue with macro arguments without rising graphic window

            _show_junctions = Boolean.valueOf(Macro.getValue(Macro.getOptions(),        "junctions", String.valueOf(_show_junctions)));
            _show_endpoints = Boolean.valueOf(Macro.getValue(Macro.getOptions(),        "endpoints", String.valueOf(_show_endpoints)));
            _enable_interactive = Boolean.valueOf(Macro.getValue(Macro.getOptions(),    "interactive", String.valueOf(_enable_interactive)));
            _calculate_directions = Boolean.valueOf(Macro.getValue(Macro.getOptions(),  "calculate_directions", String.valueOf(_calculate_directions)));
            _save_midresults = Boolean.valueOf(Macro.getValue(Macro.getOptions(),       "save_midresults", String.valueOf(_save_midresults)));
            _Dlist = Macro.getValue(Macro.getOptions(),                                 "dlist", _Dlist);
//            _sigma_ratio = Float.valueOf(Macro.getValue(Macro.getOptions(),             "sigma_ratio", String.valueOf(1f/6f)));
//            _s = Float.valueOf(Macro.getValue(Macro.getOptions(),                       "s", String.valueOf(1.2)));

            _ncc_high           = Float.valueOf(Macro.getValue(Macro.getOptions(),      "ncc_high", String.valueOf(_ncc_high)));
            _ncc_low            = Float.valueOf(Macro.getValue(Macro.getOptions(),      "ncc_low", String.valueOf(_ncc_low)));
            _likelihood_high    = Float.valueOf(Macro.getValue(Macro.getOptions(),      "likelihood_high", String.valueOf(_likelihood_high)));
            _likelihood_low     = Float.valueOf(Macro.getValue(Macro.getOptions(),      "likelihood_low", String.valueOf(_likelihood_low)));
            _smoothness_high    = Float.valueOf(Macro.getValue(Macro.getOptions(),      "smoothness_high", String.valueOf(_smoothness_high)));
            _smoothness_low     = Float.valueOf(Macro.getValue(Macro.getOptions(),      "smoothness_low", String.valueOf(_smoothness_low)));
            _dsens_list         = Macro.getValue(Macro.getOptions(),                    "dsens_list", _dsens_list);

        }

        String[] dd = _Dlist.split(","); if (dd.length==0) return;
        float[] _D = new float[dd.length];
        for (int i=0; i<dd.length; i++) _D[i] = Float.valueOf(dd[i]);

        dd = _dsens_list.split(",");  if (dd.length==0) return;
        float[] _dsens = new float[dd.length];
        for (int i = 0; i < dd.length; i++) _dsens[i] = Float.valueOf(dd[i]);

        // detection
//        IJ.log("Detector... ");
        det2d = new Detector(
                ip_load,
                1.2f, // s=1.2    (hardcoded)
                _D,
                1f/6, // sigma_ratio (hardcoded)
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
        det2d.do_directions = _calculate_directions; // f
        det2d.do_endpoints = _show_endpoints;        // t
        det2d.do_junctions = _show_junctions;        // t

        det2d.run();
        Overlay[] detection_overlays = det2d.saveDetection();  // export .det with detections

            ImageStack out_stack = new ImageStack(ip_load.getWidth(), ip_load.getHeight());

            Overlay    out_ovl = new Overlay(); // allocate joint Overlay for the stack

            for (int i = 0; i < detection_overlays.length; i++) {

                if (detection_overlays[i].size()>0) {

                    out_stack.addSlice("DET,dsens="+IJ.d2s(_dsens[i],2), ip_load.getProcessor().duplicate());

                    for (int j = 0; j < detection_overlays[i].size(); j++) {

                        Roi get_roi = detection_overlays[i].get(j);
                        get_roi.setPosition(out_stack.getSize());
                        out_ovl.add(get_roi);

                    }

                }
                else IJ.log("Empty overlay with detections.");

            }

        if (_enable_interactive) {
            ip_load.show();
            ip_load.setTitle("move-click mouse");

            // position it
            ImageWindow ip_load_iw;
            ip_load_iw = ip_load.getWindow();
            ip_load_iw.setLocation(upper_left_x+plot_w+20, upper_left_y); //+upper_left_y+plot_h+20+plot_h+20 , plot_w, 2*plot_h
            cnv = ip_load.getCanvas();
            cnv.addMouseListener(this);
            cnv.addMouseMotionListener(this);
            IJ.setTool("hand");
        }

        ImagePlus out_imp = new ImagePlus("detections", out_stack);
        out_imp.setOverlay(out_ovl);
        out_imp.show();
        ImageWindow out_imp_iw = out_imp.getWindow();
        out_imp_iw.setLocation(upper_left_x+plot_w+20 + 50, upper_left_y + 50); // plot_w, 2*plot_h

        IJ.log("done.");

    }

    public void mouseClicked(MouseEvent e)
    {
//        mouseMoved(e);
        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

/*        pfl_is3 = Delineator.getClusterMask(clickX, clickY);
        if (pfl_is3.getSize()>0) {
            pfl_im3.setStack("cluster_mask", pfl_is3);

            if (first3) {
                pfl_im3.show();
                pfl_iw3 = pfl_im3.getWindow();
                pfl_iw3.setLocation(upper_left_x, upper_left_y+upper_left_y+plot_h+50+plot_h+50);
                first3 = false;
            }
            pfl_im3.updateAndDraw();
        }*/

        pfl_is2 = (det2d!=null)? det2d.getFuzzyLogicResponse(clickX, clickY) : null;
        if (pfl_is2.getSize()>0) { //pfl_is2!=null &&
            pfl_im2.setStack("fuzzy_logic", pfl_is2);

            if (first2) {
                pfl_im2.show();
                pfl_iw2 = pfl_im2.getWindow();
                pfl_iw2.setLocation(upper_left_x, upper_left_y+plot_h+50);
                first2 = false;
            }

//            pfl_im2.setStack("fuzzy_logic", pfl_is2);
            pfl_im2.updateAndDraw();
        }

        IJ.setTool("hand");

    }

    public void mouseMoved(MouseEvent e) {
        //mouseClicked(e);

        int clickX = cnv.offScreenX(e.getX());
        int clickY = cnv.offScreenY(e.getY());

        Overlay ov_to_add = Delineator.getDelineationOverlay(clickX, clickY);
        cnv.setOverlay(ov_to_add);

        Detector.print(clickX, clickY); // print extracted features
        pfl_is4 = PeakExtractor.getProfileWithPeaks(clickX, clickY);

        if (pfl_is4.getSize()>0)
            pfl_im4.setStack("angular_profile", pfl_is4);

        if (first4) {
            pfl_im4.show();
            pfl_iw4 = pfl_im4.getWindow();
            pfl_iw4.setLocation(upper_left_x, upper_left_y);
            first4 = false;
        }
        pfl_im4.updateAndDraw();

        IJ.setTool("hand");

    }

    public void mousePressed(MouseEvent e)  {}
    public void mouseReleased(MouseEvent e) {}
    public void mouseEntered(MouseEvent e)  {}
    public void mouseExited(MouseEvent e)   {}
    public void mouseDragged(MouseEvent e)  {}

}
