package com.braincadet.npinpoint;

import ij.ImagePlus;
import ij.io.FileSaver;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;

import java.io.*;
import java.util.Arrays;

/**
 * Threaded implementation of the masker module - robust foreground point extraction
 * Separates foreground from the background to reduce the computation by comparing and thresholding the difference between
 * 95th percentile and 5th percentile of the values within a radius surrounding some location
 * where 5th percentile is considered the background estimate and 95th percentile foreground estimate
 * added option to expand the obtained mask by using structural element
 */
public class Masker extends Thread {

    private int begN, endN;

    public static int 				image_width;
    public static int 				image_height;

    public static float[][]			inimg_xy;

    private static float            radiusCheck;
    private static float 			globalTh;
    private static int              marginPix;
    private static float            percentile;

    /*
    OUTPUT
     */
    public static byte[][]			back_xy; 	// background estimate
    public static boolean[][]		mask_xy;    // output mask
    public static float[][]  		criteria;	//
    public static int[][] 			i2xy;       // mapping
    public static  int[][]			xy2i;

    public Masker(int n0, int n1) {
        begN = n0;
        endN = n1;
    }

    public static void loadTemplate(float[][] _inimg_xy, int _margin, float _radiusNbhoodCheck, float _percentile)
    {

        radiusCheck     = _radiusNbhoodCheck;
        percentile      = (_percentile<0)? 0 : (_percentile>100)? 100 : _percentile;
        percentile      = percentile / 5;

        marginPix       = (int) Math.ceil(radiusCheck);
        marginPix = (_margin>marginPix)? _margin : marginPix ;  // narrow down selection in XY plane

        inimg_xy = _inimg_xy;
        image_height 	= inimg_xy[0].length;
        image_width 	= inimg_xy.length;

		/*
			initialize outputs
		 */
        back_xy = new byte[image_width][image_height];
        mask_xy = new boolean[image_width][image_height];

        criteria = new float[image_width][image_height];

        i2xy 	= new int[1][1];  // values known after run()
        xy2i 	= new int[1][1];  // values known after run()

    }

    public void run()
    {

        int circNeighSize = sizeCircularNbhood(radiusCheck);
        float[] circNeigh = new float[circNeighSize];

        for (int locIdx=begN; locIdx<endN; locIdx++) {

            int atX = locIdx%image_width;
            int atY = locIdx/image_width;

            boolean processIt =
                    (atX>=marginPix) && (atY>=marginPix) && //(atZ>=marginLay) &&
                            (atX<image_width-marginPix-1) && (atY<image_height-marginPix-1);// && (atZ<image_length-marginLay-1);

            if (processIt) {

                extractCircularNbhood(atX, atY, radiusCheck, circNeigh, inimg_xy); // extract from the image
                float m05 	= Stat.quantile(circNeigh, 1 , 20);
                float m95 	= Stat.quantile(circNeigh, 19, 20);

                back_xy[atX][atY] 	= (byte) Math.round(m05);
                criteria[atX][atY] 	= m95 - m05;

            }

        }

    }

    public static void defineThreshold()
    {

        float[][] criteria_local_max = new float[criteria.length][criteria[0].length];

        int circNeighSize = sizeCircularNbhood(radiusCheck);
        float[] circNeigh = new float[circNeighSize];

        for (int x=0; x<criteria.length; x++) {
            for (int y=0; y<criteria[0].length; y++) {
                extractCircularNbhood(x, y, radiusCheck, circNeigh, criteria);
                criteria_local_max[x][y] = Stat.get_max(circNeigh);
            }
        }

        for (int ii=0; ii<criteria.length; ii++) {
            for (int jj=0; jj<criteria[0].length; jj++) {
                criteria[ii][jj] = criteria_local_max[ii][jj];
            }
        }

        float[] criteria_temp = new float[criteria.length*criteria[0].length];
        int		cnt = 0;
        for (int ii=0; ii<criteria.length; ii++) {
            for (int jj=0; jj<criteria[0].length; jj++) {
                criteria_temp[cnt] = criteria_local_max[ii][jj];
                cnt++;
            }
        }

        globalTh = Stat.quantile(criteria_temp, (int) percentile, 20);

//        cnt = 0;
        for (int xx=0; xx<image_width; xx++) {
            for (int yy=0; yy<image_height; yy++) {

                if (criteria[xx][yy] > globalTh) {
                    mask_xy[xx][yy] = true;
                }
                else {
                    mask_xy[xx][yy] = false;
                }

//                criteria[xx][yy] = criteria_temp[cnt];

//                cnt++;

//				if (criteria[xx][yy] > globalTh) {
//					mask_xy[xx][yy] = true;
//				}
//				else {
//					mask_xy[xx][yy] = false;
//				}
            }
        }

    }

    public static void formRemainingOutputs()
    {

        xy2i 	= new int[image_width][image_height];

        int cnt = 0;
        //for (int zz=0; zz<image_length; zz++) {
        for (int xx=0; xx<image_width; xx++) {
            for (int yy=0; yy<image_height; yy++) {

                if (mask_xy[xx][yy]) {
                    xy2i[xx][yy] = cnt;
                    cnt++;
                }
                else {
                    xy2i[xx][yy] = -1;
                }

            }
        }
        //}

//		float perc = (cnt*100f) / (image_width*image_height);
//		IJ.log(String.format("%3.2f %% vol. foreground extracted", perc));

        // component: foreground point list: locations and background estimates
        i2xy = new int[cnt][2];
//		mo.backgroundEst = new byte[cnt];

        cnt =0;
//		for (int zz=0; zz<image_length; zz++) {
        for (int xx=0; xx<image_width; xx++) {
            for (int yy=0; yy<image_height; yy++) {

                if (mask_xy[xx][yy]) {

//						mo.foregroundLocsZXY[cnt][0] = zz;
                    i2xy[cnt][0] = xx;
                    i2xy[cnt][1] = yy;

//						mo.backgroundEst[cnt] = Masker3D.back3[zz][xx][yy];

                    cnt++;

                }

            }
        }
//		}


    }

    public static ByteProcessor getMask()
    {

        byte[] out = new byte[image_height*image_width];
        for (int i=0; i<out.length; i++) {
            out[i] = mask_xy[i%image_width][i/image_width]? (byte)255 : (byte)0;
        }
        return new ByteProcessor(image_width, image_height, out);

    }

    public static FloatProcessor getCriteria()
    {
        return new FloatProcessor(criteria);
    }

    public static ByteProcessor getBackground()
    {

        byte[] out = new byte[image_height*image_width];
        for (int i=0; i<out.length; i++) {
            out[i] = back_xy[i%image_width][i/image_width];
        }
        return new ByteProcessor(image_width, image_height, out);

    }

    private static int sizeCircularNbhood(float sphereRadius)
    {

        //int rLayers 	= Math.round(sphereRadius / zDist);
        int rPix 		= Math.round(sphereRadius);

        int cnt = 0;

        for (int a=-rPix; a<=rPix; a++) {
            for (int b=-rPix; b<=rPix; b++) {
                //for (int cLayers=-rLayers; cLayers<=rLayers; cLayers++) { // loop in layers
                //float c = cLayers * zDist;   // back to pixels
                if ( a*a+b*b <= rPix*rPix ) cnt++; // +c*c
                //}
            }
        }

        return cnt;
    }

    private static void extractCircularNbhood(int atX, int atY, float sphereRadius, float[] values, float[][] _inimg_xy)
    {

        int rPix = Math.round(sphereRadius);
        //int rLay = Math.round(sphereRadius / zDist);

        if (atX-rPix>=0 && atY-rPix>=0 && atX+rPix<_inimg_xy.length && atY+rPix<_inimg_xy[0].length) {  // && atZ-rLay>=0  // && atZ+rLay < image_length

            int cnt = 0;

            for (int xLoc=atX-rPix; xLoc<=atX+rPix; xLoc++) {
                for (int yLoc=atY-rPix; yLoc<=atY+rPix; yLoc++) {
                    //for (int zLoc=atZ-rLay; zLoc<=atZ+rLay; zLoc++) { // loop in layers
                    //float c = (zLoc-atZ) * zDist;   // back to pixels
                    if ( (xLoc-atX)*(xLoc-atX)+(yLoc-atY)*(yLoc-atY) <= rPix*rPix ) { // +c*c

                        values[cnt] = _inimg_xy[xLoc][yLoc];// instack3[atZ][xLoc][yLoc];
                        cnt++;

                    }
                    //}
                }
            }
        }
        else {
            Arrays.fill(values, 0f); // set all values to zero
        }

    }

    private static float medianAtPoint(int x, int y, float[][] _inimg_xy)
    {

        float[] nhood = new float[9];

        if (x>0 && x<_inimg_xy.length-1 && y>0 && y<_inimg_xy[0].length-1) {

            int cnt = 0;
            for (int dx=-1; dx<=1; dx++) {
                for (int dy=-1; dy<=1; dy++) {
                    nhood[cnt] = _inimg_xy[x+dx][y+dy];
                    cnt++;
                }
            }

        }

        return Stat.median(nhood);

    }

	/*
		EXPORT FILE
	 */

    public static void exportI2xyCsv(String file_path) {

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        for (int ii=0; ii<i2xy.length; ii++) {  //main loop
            for (int jj=0; jj<i2xy[ii].length; jj++) {
                logWriter.print(String.format("%6d", i2xy[ii][jj]));
                if (jj<i2xy[ii].length-1) {
                    logWriter.print(",\t");
                }
                else {
                    logWriter.print("\n");
                }
            }
        }

        logWriter.close(); // close log

    }

    public static void exportXy2iCsv(String file_path) {

        PrintWriter logWriter = null; //initialize writer

        try {
            logWriter = new PrintWriter(file_path);
            logWriter.print("");
            logWriter.close();
        } catch (FileNotFoundException ex) {}   // empty the file before logging...

        try {                                   // initialize detection log file
            logWriter = new PrintWriter(new BufferedWriter(new FileWriter(file_path, true)));
        } catch (IOException e) {}

        for (int ii=0; ii<xy2i.length; ii++) {  //main loop
            for (int jj=0; jj<xy2i[ii].length; jj++) {
                logWriter.print(String.format("%15d", xy2i[ii][jj]));
                if (jj<xy2i[ii].length-1) {
                    logWriter.print(",\t");
                }
                else {
                    logWriter.print("\n");
                }
            }
        }

        logWriter.close(); // close log

    }

    public static void exportForegroundMask(String file_path) {

        ByteProcessor bp = getMask();
        ImagePlus ip = new ImagePlus("foreground_mask", bp);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(file_path);

    }

    public static void exportEstBackground(String file_path) {

        ByteProcessor bp = getBackground();
        ImagePlus ip = new ImagePlus("est_background", bp);
        FileSaver fs = new FileSaver(ip);
        fs.saveAsTiff(file_path);

    }

}
